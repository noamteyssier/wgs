#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.spatial.distance import squareform

sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})


def tri(df):
    ndf = pd.DataFrame(np.tril(df))
    ndf.index = df.index
    ndf.columns = df.columns
    return ndf

def main():
    snp_fn = "data/snpDiff.tab"
    snps = pd.read_csv(snp_fn, sep="\t")

    inds = snps[snps.sample_a == snps.sample_b].loc[:,['shared','sample_a']]
    sns.distplot(inds.shared, kde=False)


    piv = pd.pivot_table(
        snps,
        index='sample_a',
        columns='sample_b',
        values=['uniq_a', 'shared']
    )


    shared = tri(piv.loc[:, 'shared'])
    unique = tri(piv.loc[:, 'uniq_a'])

    shared
    unique

    sns.heatmap(unique)

    pass


if __name__ == '__main__':
    main()
