# The script is used to integrate the boxplot of each modification on the pearson correlation coefficient

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def integrate_boxplot(mod, path_global, path_global_all, path_global_pair, path_outTSS_all, path_outTSS_pair, path_without_h3k4me3):
    df_global = pd.read_csv(path_global, header=0)
    df_global_all = pd.read_csv(path_global_all, header=0)
    df_global_pair = pd.read_csv(path_global_pair, header=0)
    df_outTSS_all = pd.read_csv(path_outTSS_all, header=0)
    df_outTSS_pair = pd.read_csv(path_outTSS_pair, header=0)
    df_without_h3k4me3 = pd.read_csv(path_without_h3k4me3, header=0)
    ndf_global_all = pd.DataFrame(np.array(df_global_all['drm_pcc']), index=[df_global_all['chrom'], df_global_all['start']], columns=['global_all'])
    ndf_global_pair = pd.DataFrame(np.array(df_global_pair['drm_pcc']), index=[df_global_pair['chrom'], df_global_pair['start']], columns=['global_pair'])
    ndf_outTSS_all = pd.DataFrame(np.array(df_outTSS_all['drm_pcc']), index=[df_outTSS_all['chrom'], df_outTSS_all['start']], columns=['outTSS_all'])
    ndf_outTSS_pair = pd.DataFrame(np.array(df_outTSS_pair['drm_pcc']), index=[df_outTSS_pair['chrom'], df_outTSS_pair['start']], columns=['outTSS_pair'])
    ndf_without_h3k4me3 = pd.DataFrame(np.array(df_without_h3k4me3['drm_pcc']), index=[df_without_h3k4me3['chrom'], df_without_h3k4me3['start']], columns=['out_H3K4me3'])
    df_integrate = pd.DataFrame(index=[df_global['chrom'], df_global['start']])
    column_list = ['global_all', 'global_pair', 'outTSS_all', 'outTSS_pair', 'out_H3K4me3']
    ndf_list = [ndf_global_all, ndf_global_pair, ndf_outTSS_all, ndf_outTSS_pair, ndf_without_h3k4me3]
    for column, ndf in zip(column_list, ndf_list):
        df_integrate[column] = ndf[column]
    #df_integrate.to_csv('D:\\G_project\\Zhen_G\\data\\CorrelationC\\integrated\\{0}.csv'.format(mod))
    df_integrate.boxplot(column=column_list)
    plt.title(mod)
    plt.show()

Mod_list = ['H3K4me1_H3K4me3_drm_PCC', 'H3K4me1_H3K27ac_drm_PCC', 'H3K4me3_H3K27ac_drm_PCC']
Path_global = 'D:\\G_project\\Zhen_G\\data\\CorrelationC\\global\\drm_pcc.csv'
Path_global_all_list = ['D:\\G_project\\Zhen_G\\data\\CorrelationC\\global\\peak_all\\H3K4me1_H3K4me3_CC.csv',
                        'D:\\G_project\\Zhen_G\\data\\CorrelationC\\global\\peak_all\\H3K4me1_H3K27ac_CC.csv',
                        'D:\\G_project\\Zhen_G\\data\\CorrelationC\\global\\peak_all\\H3K4me3_H3K27ac_CC.csv']
Path_global_pair_list = ['D:\\G_project\\Zhen_G\\data\\CorrelationC\\global\\peak_and\\H3K4me1_H3K4me3_CC.csv',
                         'D:\\G_project\\Zhen_G\\data\\CorrelationC\\global\\peak_and\\H3K4me1_H3K27ac_CC.csv',
                         'D:\\G_project\\Zhen_G\\data\\CorrelationC\\global\\peak_and\\H3K4me3_H3K27ac_CC.csv']
Path_outTSS_all_list = ['D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me1_H3K4me3_CC.csv',
                        'D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me1_H3K27ac_CC.csv',
                        'D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_all\\enhancer_H3K4me3_H3K27ac_CC.csv']
Path_outTSS_pair_list = ['D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_pair\\enhancer_H3K4me1_H3K4me3_CC.csv',
                         'D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_pair\\enhancer_H3K4me1_H3K27ac_CC.csv',
                         'D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_peak_pair\\enhancer_H3K4me3_H3K27ac_CC.csv']
Path_without_h3k4me3_list = ['D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_without_h3k4me3\\enhancer_H3K4me1_H3K4me3_CC.csv',
                             'D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_without_h3k4me3\\enhancer_H3K4me1_H3K27ac_CC.csv',
                             'D:\\G_project\\Zhen_G\\data\\CorrelationC\\enhancer\\enhancer_region_without_h3k4me3\\enhancer_H3K4me3_H3K27ac_CC.csv']

for Mod, Path_global_all, Path_global_pair, Path_outTSS_all, Path_outTSS_pair, Path_without_h3k4me3 in zip(Mod_list, Path_global_all_list, Path_global_pair_list, Path_outTSS_all_list, Path_outTSS_pair_list, Path_without_h3k4me3_list):
    integrate_boxplot(Mod, Path_global, Path_global_all, Path_global_pair, Path_outTSS_all, Path_outTSS_pair, Path_without_h3k4me3)


