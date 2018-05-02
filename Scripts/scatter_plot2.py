# The script is used to draw a scatter plot based on the pearson correlation coefficients between three histone modifications and gene expression

import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def correlation_cal(raw_file1, raw_file2):
    df_file1 = pd.read_table(raw_file1, sep='\t')
    df_file2 = pd.read_table(raw_file2, sep='\t')
    df_correlation = df_file1.loc[:, ['chrom', 'start', 'end']]
    if len(df_file1) != len(df_file2):
        print("File length not equal!")
    else:
        for i in range(len(df_file1)):
            pcc_tup = stats.pearsonr(df_file1.iloc[i, 3:49], df_file2.iloc[i, 3:49])
            scc_tup = stats.spearmanr(df_file1.iloc[i, 3:49], df_file2.iloc[i, 3:49])
            df_correlation.loc[i, 'pcc'] = pcc_tup[0]
            df_correlation.loc[i, 'pcc_Pvalue'] = pcc_tup[1]
            df_correlation.loc[i, 'scc'] = scc_tup[0]
            df_correlation.loc[i, 'scc_Pvalue'] = scc_tup[1]
    return df_correlation
'''
H3K4me1_file = r'D:\G_project\Zhen_G\data\LCLs_C2.C1Y12_ref_norm_signals\46samples\LCLs_C2.H3K4me1_drm_norm_signals.txt'
H3K4me3_file = r'D:\G_project\Zhen_G\data\LCLs_C2.C1Y12_ref_norm_signals\46samples\LCLs_C2.H3K4me3_drm_norm_signals.txt'
H3K27ac_file = r'D:\G_project\Zhen_G\data\LCLs_C2.C1Y12_ref_norm_signals\46samples\LCLs_C2.H3K27ac_drm_norm_signals.txt'
df_H3K4me1_H3K4me3_cc = correlation_cal(H3K4me1_file, H3K4me3_file)
df_H3K4me1_H3K27ac_cc = correlation_cal(H3K4me1_file, H3K27ac_file)
df_H3K4me3_H3K27ac_cc = correlation_cal(H3K4me3_file, H3K27ac_file)
df_H3K4me1_H3K4me3_cc.to_csv(r'D:\G_project\Zhen_G\data\CorrelationC\46samples_cc\H3K4me1_H3K4me3_cc.txt', sep='\t', index=None)
df_H3K4me1_H3K27ac_cc.to_csv(r'D:\G_project\Zhen_G\data\CorrelationC\46samples_cc\H3K4me1_H3K27ac_cc.txt', sep='\t', index=None)
df_H3K4me3_H3K27ac_cc.to_csv(r'D:\G_project\Zhen_G\data\CorrelationC\46samples_cc\H3K4me3_H3K27ac_cc.txt', sep='\t', index=None)
print("Calcualte histone modification cc finished!")


df_H3K4me1_gene_cc = pd.read_table(r'D:\G_project\Zhen_G\data\expression_correlation\CC\with_H3K4me3\H3K4me1_p2g_cc.txt')
df_H3K4me3_gene_cc = pd.read_table(r'D:\G_project\Zhen_G\data\expression_correlation\CC\with_H3K4me3\H3K4me3_p2g_cc.txt')
df_H3K27ac_gene_cc = pd.read_table(r'D:\G_project\Zhen_G\data\expression_correlation\CC\with_H3K4me3\H3K27ac_p2g_cc.txt')
df_integrate = df_H3K4me1_gene_cc.loc[:, ['peak_chrom', 'peak_start', 'peak_end', 'gene_name', 'gene_name2', 'TSS_site']]
df_integrate.index = [df_integrate['peak_chrom'], df_integrate['peak_start']]
for Mod, df_cc in zip(['H3K4me1_H3K4me3', 'H3K4me1_H3K27ac', 'H3K4me3_H3K27ac'], [df_H3K4me1_H3K4me3_cc, df_H3K4me1_H3K27ac_cc, df_H3K4me3_H3K27ac_cc]):
    df_cc.index = [df_cc['chrom'], df_cc['start']]
    df_integrate['{0}_pcc'.format(Mod)] = df_cc['pcc']
    df_integrate['{0}_pcc_Pvalue'.format(Mod)] = df_cc['pcc_Pvalue']
    df_integrate['{0}_scc'.format(Mod)] = df_cc['scc']
    df_integrate['{0}_scc_Pvalue'.format(Mod)] = df_cc['scc_Pvalue']
for Mod2, df_gcc in zip(['H3K4me1', 'H3K4me3', 'H3K27ac'], [df_H3K4me1_gene_cc, df_H3K4me3_gene_cc, df_H3K27ac_gene_cc]):
    df_gcc.index = [df_gcc['peak_chrom'], df_gcc['peak_start']]
    df_integrate['{0}_gene_pcc'.format(Mod2)] = df_gcc['{0}_pcc'.format(Mod2)]
    df_integrate['{0}_gene_pcc_Pvalue'.format(Mod2)] = df_gcc['{0}_pcc_Pvalue'.format(Mod2)]
    df_integrate['{0}_gene_scc'.format(Mod2)] = df_gcc['{0}_scc'.format(Mod2)]
    df_integrate['{0}_gene_scc_Pvalue'.format(Mod2)] = df_gcc['{0}_scc_Pvalue'.format(Mod2)]
df_integrate.index = range(len(df_integrate))
df_integrate.to_csv(r'D:\G_project\Zhen_G\data\expression_correlation\integrated_all.txt', sep='\t', index=None)
'''

df_integrate = pd.read_csv(r'D:\G_project\Zhen_G\data\expression_correlation\with_H3K4me3_integrated_all.txt', sep='\t')
df_integrate.plot(kind='scatter', x='H3K4me1_H3K4me3_pcc', y='H3K4me3_H3K27ac_pcc', s=5, c='H3K4me3_gene_pcc', cmap=plt.cm.seismic)
print(stats.pearsonr(df_integrate['H3K4me1_H3K4me3_pcc'], df_integrate['H3K4me3_H3K27ac_pcc']))
plt.title('H3K4me3 PCC with H3K4me1, H3K27ac and gene expression')
plt.show()

df_integrate.plot(kind='scatter', x='H3K4me1_H3K27ac_pcc', y='H3K4me3_H3K27ac_pcc', s=5, c='H3K4me3_gene_pcc', cmap=plt.cm.seismic)
print(stats.pearsonr(df_integrate['H3K4me1_H3K27ac_pcc'], df_integrate['H3K4me3_H3K27ac_pcc']))
plt.title('H3K27ac PCC with H3K4me1, H3K4me3 and H3K4me3 PCC with gene expression')
plt.show()

df_integrate.plot(kind='scatter', x='H3K4me1_H3K27ac_pcc', y='H3K4me1_H3K4me3_pcc', s=5, c='H3K4me3_gene_pcc', cmap=plt.cm.seismic)
print(stats.pearsonr(df_integrate['H3K4me1_H3K27ac_pcc'], df_integrate['H3K4me1_H3K4me3_pcc']))
plt.title('H3K4me1 PCC with H3K27ac, H3K4me3 and H3K4me3 PCC with gene expression')
plt.show()