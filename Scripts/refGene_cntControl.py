# The script is used to screen the input refGene file and generate a new refGene file for shortest distance calculation.
# The gene count control is according to the input gene counts normalization table.

import pandas as pd
import numpy as np

def refGene_cntCtrl(refGeneFile, genecnt_nmlfile):
    refGene_header = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    df_refGene = pd.read_table(refGeneFile, sep='\t', header=None, names=refGene_header)
    df_genecnt = pd.read_table(genecnt_nmlfile, sep='\t')
    df_refGene.index = df_refGene['name2']
    df_genecnt.index = df_genecnt['name2']
    df_refGene['flag'] = df_genecnt['name2']
    df_refGene = df_refGene.dropna(axis=0, how='any')
    df_refGene.index = range(len(df_refGene))
    del df_refGene['flag']
    return df_refGene

RefGeneFile = r'D:\G_project\Zhen_G\data\hg19_refGene.txt'
Genecnt_nmlfile = r'D:\G_project\Zhen_G\data\gene_counts_table_nml_withstrand.txt'
Df_refGene = refGene_cntCtrl(RefGeneFile, Genecnt_nmlfile)
Df_refGene.to_csv(r'D:\G_project\Zhen_G\data\hg19_refGene_cntControl.txt2', index=None, sep='\t', header=None)