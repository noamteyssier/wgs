#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import kurtosis, skew

sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})

def main():
    collapsed_fn = "data/collapsed_replicates.tab.gz"
    collapsed = pd.read_csv(collapsed_fn, sep='\t')


    # pivot replicates
    reps = pd.pivot_table(
            collapsed,
            index=['density', 'extraction', 'swga', 'digest', 'chrom', 'pos'],
            columns='rep',
            values='depth').\
        reset_index()


    c = reps.groupby(["density", "extraction", "swga", "digest", "chrom"]).corr().reset_index()
    depth_corr = c[c.rep == 2].loc[:,['density', 'extraction', 'swga', 'digest', 'chrom', 1]].rename(columns={1 : 'dc'})
    depth_corr = depth_corr[~((depth_corr.density == 100) & (depth_corr.extraction == 'Che') & (depth_corr.swga == 'Sof'))]


    # measure correlation
    sns.catplot(
        x='density', y = 'dc',
        kind='box', col='swga',
        row='extraction', height=10,
        hue='digest', data=depth_corr)


    skewness= collapsed.\
        groupby(['density', 'extraction', 'swga', 'digest', 'chrom', 'rep']).\
        apply(lambda x : skew(x.depth)).\
        reset_index().\
        rename(columns={0 : 'skew'})

    skewness = skewness[~((skewness.density == 100) & (skewness.extraction == 'Che') & (skewness.swga == 'Sof'))]

    # measure skewness
    sns.catplot(
        x='density', y = 'skew',
        col='swga', row='extraction',
        height=10, hue='digest',
        kind='box',
        data=skewness)


if __name__ == '__main__':
    main()
