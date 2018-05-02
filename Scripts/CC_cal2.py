# The script is used to calculate the pearsom and spearman correlation coefficient for specific file

import numpy as np
from scipy import stats

def parse_xls_to_matrix(xls_filename):
    with open(xls_filename) as file:
        matrix = []
        lines = file.readlines()
        matrix.append(lines[0].split())
        for i in range(len(lines))[1:]:
            eachlist = lines[i].split()
            matrix.append([eachlist[0]] + [int(x) for x in eachlist[1:3]] + [float(x) for x in eachlist[3:]])
    return matrix

def findpeak(list):
    if sum(list) == 0:
        return [0]
    else:
        return [1]

# m_H3K4me1_drm = parse_xls_to_matrix('..\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K4me1_drm_norm_signals.xls')
# m_H3K4me3_drm = parse_xls_to_matrix('..\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K4me3_drm_norm_signals.xls')
# m_H3K27ac_drm = parse_xls_to_matrix('..\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K27ac_drm_norm_signals.xls')
m_H3K4me1_ndrm = parse_xls_to_matrix('..\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K4me1_ndrm_norm_signals.xls')
m_H3K4me3_ndrm = parse_xls_to_matrix('..\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K4me3_ndrm_norm_signals.xls')
m_H3K27ac_ndrm = parse_xls_to_matrix('..\data\LCLs_C2.C1Y12_ref_norm_signals\LCLs_C2.H3K27ac_ndrm_norm_signals.xls')

with open('ndrm_scc.xls', 'w') as file:
    lines = ['chrom\tstart\tend\tH3K4me1_H3K4me3_SCC\tH3K4me1_H3K4me3_SCC_P_value\tH3K4me1_H3K27ac_SCC\tH3K4me1_H3K27ac_SCC_P_value\tH3K4me3_H3K27ac_SCC\tH3K4me3_H3K27ac_SCC_P_value\tH3K4me1_peak\tH3K4me3_peak\tH3K27ac_peak\n']
    for i in range(len(m_H3K4me1_ndrm))[1:]:
        H3K4me1_peak = findpeak(m_H3K4me1_ndrm[i][50:97])
        H3K4me3_peak = findpeak(m_H3K4me3_ndrm[i][50:97])
        H3K27ac_peak = findpeak(m_H3K27ac_ndrm[i][50:97])
        line = m_H3K4me1_ndrm[i][0:3] + list(stats.spearmanr(m_H3K4me1_ndrm[i][3:50], m_H3K4me3_ndrm[i][3:50])) + list(stats.spearmanr(m_H3K4me1_ndrm[i][3:50], m_H3K27ac_ndrm[i][3:50])) + list(stats.spearmanr(m_H3K4me3_ndrm[i][3:50], m_H3K27ac_ndrm[i][3:50])) + H3K4me1_peak + H3K4me3_peak + H3K27ac_peak
        lines.append('\t'.join([str(x) for x in line]) + '\n')
    file.writelines(lines)