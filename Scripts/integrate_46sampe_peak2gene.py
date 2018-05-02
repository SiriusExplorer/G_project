# The script is used to integrate the 46 people's normalized histone modification peaks under specific condition and the corresponed gene.
# The script will calculate the Pearson and the Spearman correlation between specific histone modification and gene expression

import pandas as pd
import numpy as np
from scipy import stats

def find_peak(list):
    if sum(list) == 0:
        return 0
    else:
        return 1

def sample46_hispeak(h3k4me1_file, h3k4me3_file, h3k27ac_file, gene_count_file):
    df_h3k4me1 = pd.read_csv(h3k4me1_file, sep='\t')
    df_h3k4me3 = pd.read_csv(h3k4me3_file, sep='\t')
    df_h3k27ac = pd.read_csv(h3k27ac_file, sep='\t')
    df_genecount = pd.read_csv(gene_count_file, sep='\t')
    df_integrate_peaks = df_h3k4me1.loc[:, ['chrom', 'start', 'end']]
    genecount_sample_list = df_genecount.columns[1:47]
    h3k4me1_peak_name_list = [x.replace('_gene_cnt_nml', '_{0}.occupancy'.format('H3K4me1')) for x in genecount_sample_list]
    h3k4me3_peak_name_list = [x.replace('_gene_cnt_nml', '_{0}.occupancy'.format('H3K4me3')) for x in genecount_sample_list]
    h3k27ac_peak_name_list = [x.replace('_gene_cnt_nml', '_{0}.occupancy'.format('H3K27ac')) for x in genecount_sample_list]
    for i in range(len(df_integrate_peaks)):
        df_integrate_peaks.loc[i, 'H3K4me1_peak'] = find_peak(df_h3k4me1.loc[i, h3k4me1_peak_name_list])
        df_integrate_peaks.loc[i, 'H3K4me3_peak'] = find_peak(df_h3k4me3.loc[i, h3k4me3_peak_name_list])
        df_integrate_peaks.loc[i, 'H3K27ac_peak'] = find_peak(df_h3k27ac.loc[i, h3k27ac_peak_name_list])
    return df_integrate_peaks

def integrate(mod, peak2gene_file, histone_sample_file, gene_count_file, df_integrate_peaks):
    df_peak2gene = pd.read_csv(peak2gene_file, sep='\t')
    df_hisample = pd.read_csv(histone_sample_file, sep='\t')
    df_genecount = pd.read_csv(gene_count_file, sep='\t')
    genecount_sample_list = list(df_genecount.columns[1:47])
    hisample_name_list = [x.replace('_gene_cnt_nml', '_{0}.read_cnt'.format(mod)) for x in genecount_sample_list]
    df_integrate = df_peak2gene.loc[:, ['peak_chrom', 'peak_start', 'peak_end', 'gene_name', 'gene_name2', 'TSS_site']]
    df_integrate.index = [df_integrate['peak_chrom'], df_integrate['peak_start']]
    df_integrate_peaks.index = [df_integrate_peaks['chrom'], df_integrate_peaks['start']]
    df_peak2gene.index = [df_peak2gene['peak_chrom'], df_peak2gene['peak_start']]
    df_hisample.index = [df_hisample['chrom'], df_hisample['start']]
    df_genecount.index = df_genecount['name2']
    for sample_name in hisample_name_list:
        df_integrate[sample_name] = df_hisample[sample_name]
    df_integrate['H3K4me1_peak'] = df_integrate_peaks['H3K4me1_peak']
    df_integrate['H3K4me3_peak'] = df_integrate_peaks['H3K4me3_peak']
    df_integrate['H3K27ac_peak'] = df_integrate_peaks['H3K27ac_peak']
    df_integrate = df_integrate[(df_integrate['H3K4me1_peak'] == 1) & (df_integrate['H3K4me3_peak'] != 1) & (df_integrate['H3K27ac_peak'] == 1)]
    df_integrate.index = [df_integrate['gene_name2']]
    for genecount_sample in genecount_sample_list:
        df_integrate[genecount_sample] = df_genecount[genecount_sample]
    df_integrate.index = range(len(df_integrate))
    return df_integrate


H3K4me1_file = r'D:\G_project\Zhen_G\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K4me1_drm_norm_signals.xls'
H3K4me3_file = r'D:\G_project\Zhen_G\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K4me3_drm_norm_signals.xls'
H3K27ac_file = r'D:\G_project\Zhen_G\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K27ac_drm_norm_signals.xls'
Gene_count_file = r'D:\G_project\Zhen_G\data\gene_counts_table_nml_sim.txt'
Peak2gene_file = r'D:\G_project\Zhen_G\data\expression_correlation\peak_region2gene\ref_cntControl_peak2gene\out_H3K4me3_Cdistance_cntCtrl_binselected.txt'
Df_integrate_peaks = sample46_hispeak(H3K4me1_file, H3K4me3_file, H3K27ac_file, Gene_count_file)
for Mod, Histone_file in zip(['H3K4me1', 'H3K4me3', 'H3K27ac'], [H3K4me1_file, H3K4me3_file, H3K27ac_file]):
    Df_integrate = integrate(Mod, Peak2gene_file, Histone_file, Gene_count_file, Df_integrate_peaks)
    Df_integrate.to_csv(r'D:\G_project\Zhen_G\data\expression_correlation\46sample_nml_peak2gene\ref_cntControl_peak2gene\out_H3K4me3\{0}_nml_peak2genecnt.txt'.format(Mod), sep='\t', index=None)



