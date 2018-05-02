# The script is used to check the peak in drm and ndrm files

import matplotlib.pyplot as plt
import pandas as pd
from pandas import Series, DataFrame

df_drm_pcc = pd.read_csv('D:\G_project\Zhen_G\data\CorrelationC\csv\drm_pcc.csv', header = 0)
df_drm_scc = pd.read_csv('D:\G_project\Zhen_G\data\CorrelationC\csv\drm_scc.csv', header = 0)
df_ndrm_pcc = pd.read_csv('D:\G_project\Zhen_G\data\CorrelationC\csv\\ndrm_pcc.csv', header = 0)
df_ndrm_scc = pd.read_csv('D:\G_project\Zhen_G\data\CorrelationC\csv\\ndrm_scc.csv', header = 0)

print(sum(df_drm_pcc['H3K4me1_peak'] != df_ndrm_pcc['H3K4me1_peak']))
print(sum(df_drm_pcc['H3K4me3_peak'] != df_ndrm_pcc['H3K4me3_peak']))
print(sum(df_drm_pcc['H3K27ac_peak'] != df_ndrm_pcc['H3K27ac_peak']))

print(sum(df_drm_scc['H3K4me1_peak'] != df_ndrm_scc['H3K4me1_peak']))
print(sum(df_drm_scc['H3K4me3_peak'] != df_ndrm_scc['H3K4me3_peak']))
print(sum(df_drm_scc['H3K27ac_peak'] != df_ndrm_scc['H3K27ac_peak']))