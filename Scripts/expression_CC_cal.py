# The script is used to calculate the correlation coefficient between the gene counts and histone modification signal
# The input normalization data file should have 46 samples of histone modification signal and gene count (normalized)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def pcc_scc_cal(mod, nml_peak2genecnt_file):
    df_nml_p2g = pd.read_csv(nml_peak2genecnt_file, sep='\t')
    df_cc = df_nml_p2g.loc[:, ['peak_chrom', 'peak_start', 'peak_end', 'gene_name', 'gene_name2', 'TSS_site']]
    for i in range(len(df_nml_p2g)):
        pcc_tup = stats.pearsonr(df_nml_p2g.iloc[i, 6:52], df_nml_p2g.iloc[i, 55:101])
        scc_tup = stats.spearmanr(df_nml_p2g.iloc[i, 6:52], df_nml_p2g.iloc[i, 55:101])
        df_cc.loc[i, '{0}_pcc'.format(mod)] = pcc_tup[0]
        df_cc.loc[i, '{0}_pcc_Pvalue'.format(mod)] = pcc_tup[1]
        df_cc.loc[i, '{0}_scc'.format(mod)] = scc_tup[0]
        df_cc.loc[i, '{0}_scc_Pvalue'.format(mod)] = scc_tup[1]
    return df_cc

Mod_list = ['H3K4me1', 'H3K4me3', 'H3K27ac']
for Mod in Mod_list:
    Nml_peak2genecnt_file = r'D:\G_project\Zhen_G\data\expression_correlation\46sample_nml_peak2gene\ref_cntControl_peak2gene\with_H3K4me3\{0}_nml_peak2genecnt.txt'.format(Mod)
    Df_cc = pcc_scc_cal(Mod, Nml_peak2genecnt_file)
    Df_cc.to_csv(r'D:\G_project\Zhen_G\data\expression_correlation\CC\with_H3K4me3\{0}_p2g_cc.txt'.format(Mod), sep='\t', index=None)
    Df_cc.boxplot(column=['{0}_pcc'.format(Mod), '{0}_scc'.format(Mod)])
    plt.title('{0} Pearson & Spearman CC with H3K4me3'.format(Mod))
    plt.show()
