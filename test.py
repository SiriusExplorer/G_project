import xlrd
import xlwt
from scipy import stats
import numpy as np

filelines = []
table_H3K4me3 = xlrd.open_workbook('LCLs_C2.H3K4me3_ndrm_norm_signals.xls').sheets()[0]
table_H3K4me1 = xlrd.open_workbook('LCLs_C2.H3K4me1_ndrm_norm_signals.xls').sheets()[0]
filelines.append(['chrom', 'start', 'end', 'pearsonr', 'p_value'])
for i in [x for x in range(600) if x != 0]:
    part2 = list(stats.pearsonr(table_H3K4me1.row_values(i)[3:], table_H3K4me3.row_values(i)[3:]))
    part1 = list(table_H3K4me1.row_values(i)[0:3])
    filelines.append(part1 + part2)
matrix = np.array(filelines)
np.save('test.npy', matrix)
