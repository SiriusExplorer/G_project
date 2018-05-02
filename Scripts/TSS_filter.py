# The script is used to filter peak region out of TSS +- 5kb

import pandas as pd
import numpy as np

def parse_tss_region(refGeneFile): # Return a list of TSS region TSS site +- 5kb
    TSS_region_matrix = []
    refGene_header = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    df_refGene = pd.read_table(refGeneFile, sep='\t', header=None, names=refGene_header)
    for i in range(len(df_refGene)):
        if df_refGene.loc[i, 'strand'] == '+':
            TSS_site = int(df_refGene.loc[i, 'txStart'])
        elif df_refGene.loc[i, 'strand'] == '-':
            TSS_site = int(df_refGene.loc[i, 'txEnd'])
        else:
            print("RefGene file column strand with wrong format, line:{0}".format(i+1))
        TSS_region_matrix.append([TSS_site-5000+1, TSS_site+5000+1])
    return np.array(TSS_region_matrix)

def check_overlap(peak_region_line, TSS_region_matrix): # Check whether one line of peak region overlap with TSS region matrix
    for i in range(TSS_region_matrix.shape[0]):
        if TSS_region_matrix[i, 0] <= peak_region_line[0] & peak_region_line[0] <= TSS_region_matrix[i, 1]:
            return 0
        elif TSS_region_matrix[i, 0] <= peak_region_line[1] & peak_region_line[1] <= TSS_region_matrix[i, 1]:
            return 0
        elif peak_region_line[0] <= TSS_region_matrix[i, 0] & TSS_region_matrix[i, 1] <= peak_region_line[1]:
            return 0
    if i != TSS_region_matrix.shape[0] - 1:
        print("Error in the overlap checking")
    return 1

TSS_region = parse_tss_region('D:\\G_project\\Zhen_G\\data\\hg19_refGene.txt')
df_CC = pd.read_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\global\\drm_pcc.csv', header=0)
flag_list = []
peak_region = np.array(df_CC.loc[:, ['start', 'end']])
for i in range(peak_region.shape[0]):
    if check_overlap(peak_region[i], TSS_region) == 1:
        flag_list.append(True)
    else:
        flag_list.append(False)
df_CC_outTSS = df_CC[flag_list]
df_CC_outTSS.to_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\out_TSS\\out_TSS_drm_pcc.csv', index=False)





