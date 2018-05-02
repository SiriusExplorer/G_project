# The script is used to select the bins whose ends are same to another bin.
# For these kind of continuous bins, we only keep the first one as the enhancer region.

import pandas as pd
import numpy as np

def bin_select(peak2gene_file):
    df_peak2gene = pd.read_csv(peak2gene_file, sep='\t')
    df_peak2gene = df_peak2gene.sort_values(by=['peak_chrom', 'peak_start'])
    for i in range(len(df_peak2gene)):
        flag = 1
        if i == 0:
            df_peak2gene.loc[i, 'continuous_flag'] = flag
        else:
            if (df_peak2gene.loc[i-1, 'peak_end'] == df_peak2gene.loc[i, 'peak_start']) & (df_peak2gene.loc[i-1, 'peak_chrom'] == df_peak2gene.loc[i, 'peak_chrom']):
                flag = 0
            else:
                flag = 1
            df_peak2gene.loc[i, 'continuous_flag'] = flag
    df_peak2gene_selected = df_peak2gene[df_peak2gene['continuous_flag'] == 1]
    del df_peak2gene_selected['continuous_flag']
    return df_peak2gene_selected

Df_peak2gene = bin_select(r'D:\G_project\Zhen_G\data\expression_correlation\peak_region2gene\ref_cntControl_peak2gene\out_H3K4me3_region_closest_distance_cntControl.txt')
Df_peak2gene.to_csv(r'D:\G_project\Zhen_G\data\expression_correlation\peak_region2gene\ref_cntControl_peak2gene\out_H3K4me3_Cdistance_cntCtrl_binselected.txt', sep='\t', index=None)




