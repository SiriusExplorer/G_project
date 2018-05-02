# The script is used to find the genes within the enhancer site +- 500kb region

import pandas as pd
import numpy as np

def find_genes_within(refGeneFile, peak_region_file):
    refGene_header = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
                      'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    df_refgene = pd.read_table(refGeneFile, sep='\t', header=None, names=refGene_header)
    df_peak_region = pd.read_table(peak_region_file, sep='\t')
    df_1mb_region2genes = df_peak_region.loc[:, ['peak_chrom', 'peak_start', 'peak_end']]
    for i in range(df_1mb_region2genes):
        left_bound = df_1mb_region2genes.loc[i, 'peak_start'] - 500000
        right_bound = df_1mb_region2genes.loc[i, 'peak_end'] + 500000
