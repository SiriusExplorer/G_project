# The script is used to calculate the closest distance between the TSS site and the enhancer edge.
# In this case, we used the refGene file under count control, in case to avoid assign a low expression gene to the peak region.

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_tss_position(refGeneFile): # Return a dataframe of TSS of the given reference gene file
    refGene_header = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    df_refGene = pd.read_table(refGeneFile, sep='\t', header=None, names=refGene_header)
    df_TSS = df_refGene.loc[:, ['name', 'name2', 'chrom']]
    for i in range(len(df_refGene)):
        if '_' not in df_refGene.loc[i, 'chrom']:
            if df_refGene.loc[i, 'strand'] == '+':
                TSS_site = int(df_refGene.loc[i, 'txStart'])
            elif df_refGene.loc[i, 'strand'] == '-':
                TSS_site = int(df_refGene.loc[i, 'txEnd'])
            else:
                print("RefGene file column strand with wrong format, line:{0}".format(i+1))
            df_TSS.loc[i, 'TSS_site'] = TSS_site + 1
    df_TSS = df_TSS.dropna(axis=0, how='any')
    df_TSS.index = range(len(df_TSS))
    return df_TSS

def closest_distance_tss_peak_region(peak_file, refGeneFile):
    df_peak_file = pd.read_csv(peak_file, header=0)
    df_peak_region = df_peak_file.loc[:, ['chrom', 'start', 'end']]
    df_peak_region = df_peak_region.rename(columns={'chrom':'peak_chrom', 'start':'peak_start', 'end':'peak_end'})
    df_TSS_site = parse_tss_position(refGeneFile)
    df_TSS_site = df_TSS_site.rename(columns={'name':'gene_name', 'name2':'gene_name2', 'chrom':'gene_chrom'})
    df_closest_dis = df_peak_region.loc[:, ['peak_chrom', 'peak_start', 'peak_end']]
    for peak_i in range(len(df_peak_region)):
        closest_dis = 1000000000
        for tss_i in range(len(df_TSS_site)):
            if df_peak_region.loc[peak_i, 'peak_chrom'] == df_TSS_site.loc[tss_i, 'gene_chrom']:
                current_dis_left = abs(df_peak_region.loc[peak_i, 'peak_start'] - df_TSS_site.loc[tss_i, 'TSS_site'])
                current_dis_right = abs(df_peak_region.loc[peak_i, 'peak_end'] - df_TSS_site.loc[tss_i, 'TSS_site'])
                if min(current_dis_left, current_dis_right) < closest_dis:
                    closest_dis = min(current_dis_left, current_dis_right)
                    closest_name = df_TSS_site.loc[tss_i, 'gene_name']
                    closest_name2 = df_TSS_site.loc[tss_i, 'gene_name2']
                    closest_chrom = df_TSS_site.loc[tss_i, 'gene_chrom']
                    closest_TSS_site = df_TSS_site.loc[tss_i, 'TSS_site']
        df_closest_dis.loc[peak_i, 'closest_distance_peak2tss'] = closest_dis
        df_closest_dis.loc[peak_i, 'gene_name'] = closest_name
        df_closest_dis.loc[peak_i, 'gene_name2'] = closest_name2
        df_closest_dis.loc[peak_i, 'gene_chrom'] = closest_chrom
        df_closest_dis.loc[peak_i, 'TSS_site'] = closest_TSS_site
    return df_closest_dis

peak_file = 'D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me1_H3K4me3_CC.csv'
refGeneFile = 'D:\\G_project\\Zhen_G\\data\\hg19_refGene_cntControl.txt'
df_closest_dis = closest_distance_tss_peak_region(peak_file, refGeneFile)
df_closest_dis['closest_distance_peak2tss'].hist()
plt.title('shortest distance distribution with H3K4me3')
df_closest_dis.to_csv('D:\\G_project\\Zhen_G\\data\\expression_correlation\\peak_region2gene\\with_H3K4me3_region_closest_distance_cntControl.txt', sep='\t', index=None)
plt.show()
