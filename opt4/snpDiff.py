#!/usr/bin/env python3

import numpy as np
import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
#
# from scipy.spatial.distance import squareform
import vcf as pyvcf
import sys
# sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})


def tri(df):
    ndf = pd.DataFrame(np.tril(df))
    ndf.index = df.index
    ndf.columns = df.columns
    return ndf

def main():
    snp_fn = "data/snpDiff.tab"
    vcf_fn = "data/intersect_replicateVariance.vcf"
    vcf = pyvcf.Reader(open(vcf_fn, 'r'))

    vcf.samples
    vcf.infos

    for read in vcf:
        print(read)


if __name__ == '__main__':
    main()
