# The script is used to draw boxplots for H3K4me1_H3K4me3, H3K4me1_H3K27ac, H3K4me3_H3k27ac pairs in different CC and conditions
# Peaks of 3 kinds of histone modifications should be consistent in 4 CC file.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas import Series, DataFrame

df_drm_pcc = pd.read_csv('D:\G_project\Zhen_G\data\CorrelationC\global\drm_pcc.csv', header = 0)
df_drm_scc = pd.read_csv('D:\G_project\Zhen_G\data\CorrelationC\global\drm_scc.csv', header = 0)
df_ndrm_pcc = pd.read_csv('D:\G_project\Zhen_G\data\CorrelationC\global\\ndrm_pcc.csv', header = 0)
df_ndrm_scc = pd.read_csv('D:\G_project\Zhen_G\data\CorrelationC\global\\ndrm_scc.csv', header = 0)

df_H3K4me1_H3K4me3 = df_drm_pcc.loc[:, ['chrom', 'start', 'end', 'H3K4me1_peak', 'H3K4me3_peak', 'H3K27ac_peak']]
df_H3K4me1_H3K4me3.loc[:, 'drm_pcc'] = df_drm_pcc['H3K4me1_H3K4me3_PCC']
df_H3K4me1_H3K4me3.loc[:, 'drm_scc'] = df_drm_scc['H3K4me1_H3K4me3_SCC']
df_H3K4me1_H3K4me3.loc[:, 'ndrm_pcc'] = df_ndrm_pcc['H3K4me1_H3K4me3_PCC']
df_H3K4me1_H3K4me3.loc[:, 'ndrm_scc'] = df_ndrm_scc['H3K4me1_H3K4me3_SCC']
df_H3K4me1_H3K4me3_1peak = df_H3K4me1_H3K4me3[(df_H3K4me1_H3K4me3['H3K4me1_peak'] == 1) & (df_H3K4me1_H3K4me3['H3K4me3_peak'] == 1) & (df_H3K4me1_H3K4me3['H3K27ac_peak'] == 1)]
df_H3K4me1_H3K4me3_1peak.boxplot(column = ['drm_pcc', 'drm_scc', 'ndrm_pcc', 'ndrm_scc'])
plt.title('H3K4me1_H3K4me3')
plt.show()
df_H3K4me1_H3K4me3_1peak.to_csv('D:\G_project\Zhen_G\data\CorrelationC\pairs\peak_all\H3K4me1_H3K4me3_CC.csv', index = False)

df_H3K4me1_H3K27ac = df_drm_pcc.loc[:, ['chrom', 'start', 'end', 'H3K4me1_peak', 'H3K4me3_peak', 'H3K27ac_peak']]
df_H3K4me1_H3K27ac.loc[:, 'drm_pcc'] = df_drm_pcc['H3K4me1_H3K27ac_PCC']
df_H3K4me1_H3K27ac.loc[:, 'drm_scc'] = df_drm_scc['H3K4me1_H3K27ac_SCC']
df_H3K4me1_H3K27ac.loc[:, 'ndrm_pcc'] = df_ndrm_pcc['H3K4me1_H3K27ac_PCC']
df_H3K4me1_H3K27ac.loc[:, 'ndrm_scc'] = df_ndrm_scc['H3K4me1_H3K27ac_SCC']
df_H3K4me1_H3K27ac_1peak = df_H3K4me1_H3K27ac[(df_H3K4me1_H3K27ac['H3K4me1_peak'] == 1) & (df_H3K4me1_H3K27ac['H3K4me3_peak'] == 1) & (df_H3K4me1_H3K27ac['H3K27ac_peak'] == 1)]
df_H3K4me1_H3K27ac_1peak.boxplot(column = ['drm_pcc', 'drm_scc', 'ndrm_pcc', 'ndrm_scc'])
plt.title('H3K4me1_H3K27ac')
plt.show()
df_H3K4me1_H3K27ac_1peak.to_csv('D:\G_project\Zhen_G\data\CorrelationC\pairs\peak_all\H3K4me1_H3K27ac_CC.csv', index = False)

df_H3K4me3_H3K27ac = df_drm_pcc.loc[:, ['chrom', 'start', 'end', 'H3K4me1_peak', 'H3K4me3_peak', 'H3K27ac_peak']]
df_H3K4me3_H3K27ac.loc[:, 'drm_pcc'] = df_drm_pcc['H3K4me3_H3K27ac_PCC']
df_H3K4me3_H3K27ac.loc[:, 'drm_scc'] = df_drm_scc['H3K4me3_H3K27ac_SCC']
df_H3K4me3_H3K27ac.loc[:, 'ndrm_pcc'] = df_ndrm_pcc['H3K4me3_H3K27ac_PCC']
df_H3K4me3_H3K27ac.loc[:, 'ndrm_scc'] = df_ndrm_scc['H3K4me3_H3K27ac_SCC']
df_H3K4me3_H3K27ac_1peak = df_H3K4me3_H3K27ac[(df_H3K4me3_H3K27ac['H3K4me1_peak'] == 1) & (df_H3K4me3_H3K27ac['H3K4me3_peak'] == 1) & (df_H3K4me3_H3K27ac['H3K27ac_peak'] == 1)]
df_H3K4me3_H3K27ac_1peak.boxplot(column = ['drm_pcc', 'drm_scc', 'ndrm_pcc', 'ndrm_scc'])
plt.title('H3K4me3_H3K27ac')
plt.show()
df_H3K4me3_H3K27ac_1peak.to_csv('D:\G_project\Zhen_G\data\CorrelationC\pairs\peak_all\H3K4me3_H3K27ac_CC.csv', index = False)