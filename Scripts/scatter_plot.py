# The script is used to integrate the duplicates removed pearson correlation coefficient of H3K4me1-H3K27ac and H3K4me3-H3K27ac at the region out of TSS when all three peaks exist.

import pandas as pd
import matplotlib.pyplot as plt

df_H3K4me1_H3K27ac = pd.read_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me1_H3K27ac_CC.csv')
df_H3K4me1_H3K4me3 = pd.read_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me1_H3K4me3_CC.csv')
df_H3K4me3_H3K27ac = pd.read_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me3_H3K27ac_CC.csv')
df_integrate = pd.DataFrame()
df_integrate['H3K4me1_H3K4me3_PCC'] = df_H3K4me1_H3K4me3['drm_pcc']
df_integrate['H3K4me3_H3K27ac_PCC'] = df_H3K4me3_H3K27ac['drm_pcc']
df_integrate['H3K4me1_H3K27ac_PCC'] = df_H3K4me1_H3K27ac['drm_pcc']
df_integrate.plot(kind='scatter', x='H3K4me1_H3K27ac_PCC', y='H3K4me3_H3K27ac_PCC', s=5, c='H3K4me1_H3K4me3_PCC', cmap=plt.cm.seismic)
plt.title('PCC outTSS peak all H3K27ac')
plt.show()

df_integrate.plot(kind='scatter', x='H3K4me1_H3K4me3_PCC', y='H3K4me3_H3K27ac_PCC', s=5, c='H3K4me1_H3K27ac_PCC', cmap=plt.cm.seismic)
plt.title('PCC outTSS peak all H3K4me3')
plt.show()

df_integrate.plot(kind='scatter', x='H3K4me1_H3K4me3_PCC', y='H3K4me1_H3K27ac_PCC', s=5, c='H3K4me3_H3K27ac_PCC', cmap=plt.cm.seismic)
plt.title('PCC outTSS peak all H3K4me1')
plt.show()



# df_integrate_H3K27ac = pd.DataFrame()
# df_integrate_H3K27ac['H3K4me1_H3K27ac_PCC'] = df_H3K4me1_H3K27ac['drm_pcc']
# df_integrate_H3K27ac['H3K4me3_H3K27ac_PCC'] = df_H3K4me3_H3K27ac['drm_pcc']
# df_integrate_H3K27ac.plot(kind='scatter', x='H3K4me1_H3K27ac_PCC', y='H3K4me3_H3K27ac_PCC', s=5)
# plt.title('PCC_outTSS_peak_all_H3K27ac')
# plt.show()
#
# df_H3K4me1_H3K4me3 = pd.read_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me1_H3K4me3_CC.csv')
# df_H3K4me3_H3K27ac = pd.read_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me3_H3K27ac_CC.csv')
# df_integrate_H3K4me3 = pd.DataFrame()
# df_integrate_H3K4me3['H3K4me1_H3K4me3_PCC'] = df_H3K4me1_H3K4me3['drm_pcc']
# df_integrate_H3K4me3['H3K4me3_H3K27ac_PCC'] = df_H3K4me3_H3K27ac['drm_pcc']
# df_integrate_H3K4me3['H3K4me1_H3K27ac_PCC'] = df_H3K4me1_H3K27ac['drm_pcc']
# df_integrate_H3K4me3.plot(kind='scatter', x='H3K4me1_H3K4me3_PCC', y='H3K4me3_H3K27ac_PCC',
#                           s=5, c='H3K4me1_H3K27ac_PCC', cmap=plt.cm.seismic)
# plt.title('PCC_outTSS_peak_all_H3K4me3')
# plt.show()
#
# df_H3K4me1_H3K27ac = pd.read_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me1_H3K27ac_CC.csv')
# df_H3K4me1_H3K4me3 = pd.read_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me1_H3K4me3_CC.csv')
# df_integrate_H3K4me1 = pd.DataFrame()
# df_integrate_H3K4me1['H3K4me1_H3K27ac_PCC'] = df_H3K4me1_H3K27ac['drm_pcc']
# df_integrate_H3K4me1['H3K4me1_H3K4me3_PCC'] = df_H3K4me1_H3K4me3['drm_pcc']
# df_integrate_H3K4me1.plot(kind='scatter', x='H3K4me1_H3K4me3_PCC', y='H3K4me1_H3K27ac_PCC', s=5)
# plt.title('PCC_outTSS_peak_all_H3K4me1')
# plt.show()

