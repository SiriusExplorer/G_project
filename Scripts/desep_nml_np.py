# The script is used to integrate the gene count files from the output of htseq of each sample. A gene_sample counts matrix will be generated
# The normalization is according to the DESeq normalization process.

import os
import re
import pandas as pd
import numpy as np
from scipy import stats

def parse_htseq2series(hrseq_file): # Return a Series using gene name2 as index
    sample_name = re.search(r'GM_\d+', hrseq_file).group()
    column_name = sample_name.replace('_', '') + '_gene_cnt'
    df = pd.read_table(hrseq_file, sep='\t', header=None, names=['Gene', column_name])
    df_count = df[(df['Gene'] != '__no_feature') & (df['Gene'] != '__ambiguous') & (df['Gene'] != '__too_low_aQual') & (df['Gene'] != '__not_aligned') & (df['Gene'] != '__alignment_not_unique')]
    df_count = df_count.set_index('Gene')
    return pd.Series(df_count[column_name])

def parse_tss_position(refGeneFile): # Return a dataframe of TSS of the given reference gene file using gene name2 as index
    refGene_header = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    df_refGene = pd.read_table(refGeneFile, sep='\t', header=None, names=refGene_header)
    df_TSS = df_refGene.loc[:, ['name', 'name2', 'chrom', 'strand']]
    for i in range(len(df_refGene)):
        if '_' not in df_refGene.loc[i, 'chrom']:
            if df_refGene.loc[i, 'strand'] == '+':
                TSS_site = int(df_refGene.loc[i, 'txStart'])
            elif df_refGene.loc[i, 'strand'] == '-':
                TSS_site = int(df_refGene.loc[i, 'txEnd'])
            else:
                print("RefGene file column strand with wrong format, line:{0}".format(i+1))
            df_TSS.loc[i, 'TSS_site'] = TSS_site + 1
    df_TSS = df_TSS.set_index('name2')
    df_TSS = df_TSS.dropna(axis=0, how='any')
    return df_TSS

def generate_counts_matrix(hrseq_file_dir, refGeneFile): # Return a dataframe of gene counts table with numeral index
    df_count_matrix = parse_tss_position(refGeneFile)
    for file in os.listdir(hrseq_file_dir):
        if os.path.splitext(file)[1] == '.no_gene_cnts':
            file_path = os.path.join(hrseq_file_dir, file)
            sample_series = parse_htseq2series(file_path)
            df_count_matrix[sample_series.name] = sample_series
    df_count_matrix = df_count_matrix.reset_index('name2')
    return df_count_matrix

def gene_counts_control(df_count_matrix): # Remove the line which all of the gene counts less than 5; The index should be continuous from 0.
    for i in range(len(df_count_matrix)):
        flag = 1
        if sum(df_count_matrix.iloc[i, 5:51]) < 230:
            flag = 0
            for col in range(46):
                if df_count_matrix.iloc[i, col+5] >= 5:
                    flag = 1
                    break
        df_count_matrix.loc[i, 'control_flag'] = flag
    df_screened = df_count_matrix[(df_count_matrix['control_flag'] == 1)]
    df_screened.index = range(len(df_screened))
    del df_screened['control_flag']
    return df_screened

def deseq_normalization(df_cnt_matrix_control):
    df_cnt_nml = df_cnt_matrix_control.loc[:, ['name', 'name2', 'chrom', 'strand', 'TSS_site']]
    list_ref_cnt = []
    for i in range(len(df_cnt_matrix_control)):
        count_list = [int(x) for x in df_cnt_matrix_control.iloc[i, 5:51]]
        list_ref_cnt.append(stats.mstats.gmean(count_list))
    ary_ref_cnt = np.array(list_ref_cnt)
    for col in range(46):
        col_name = df_cnt_matrix_control.columns[col+5]
        col_name_nml = col_name + '_nml'
        ary_cnt = df_cnt_matrix_control.loc[:, col_name]
        ary_factor = ary_cnt / ary_ref_cnt
        factor_median = np.median(ary_factor[np.logical_not(np.isnan(ary_factor))])
        ary_cnt_nml = np.log2(ary_cnt / factor_median + 0.5)
        df_cnt_nml.loc[:, col_name_nml] = ary_cnt_nml
    return df_cnt_nml

Df_count_matrix = generate_counts_matrix('D:\\G_project\\Zhen_G\\data\\regularize_gene_counts_no_strand', 'D:\\G_project\\Zhen_G\\data\\hg19_refGene.txt')
Df_count_matrix.drop_duplicates('name2', keep='first', inplace=True)
Df_count_matrix.index = range(len(Df_count_matrix))
Df_count_matrix_counts_control = gene_counts_control(Df_count_matrix)
Df_cnt_nml = deseq_normalization(Df_count_matrix_counts_control)
Df_cnt_nml.to_csv('D:\\G_project\\Zhen_G\\data\\gene_counts_table_nml2.txt', index=None, sep='\t')
#Df_count_matrix_counts_control.to_csv('D:\\G_project\\Zhen_G\\data\\gene_counts_table_cntControl.csv', index=None)
