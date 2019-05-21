#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import kurtosis, skew

sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})

def add_sampleName(collapsed):
    """hacky method to add sampleNames"""
    collapsed['sampleName'] = (np.array([collapsed.density.astype('str'), collapsed.extraction, collapsed.digest, collapsed.swga, collapsed.rep.astype('str')]) + '-').sum(axis=0)
    collapsed['sampleName'] = collapsed['sampleName'].apply(lambda x : x.strip('-'))
    return collapsed
def pivot_replicates(collapsed):
    """bring pivot replicates to columns"""
    return pd.pivot_table(
        collapsed,
        index=['density', 'extraction', 'swga', 'digest', 'chrom', 'pos'],
        columns='rep',
        values='depth').reset_index()
def replicate_correlation(collapsed_replicates):
    """correlate replicates by group"""
    rep_corr = collapsed_replicates.\
        groupby(collapsed_replicates.columns.values[:-3].tolist())[[1,2]].\
        corr().\
        iloc[0::2][2].\
        reset_index().\
        drop(columns='rep').\
        rename(columns={2:'correlation'})
    return rep_corr
def replicate_skew(collapsed):
    """calculate skew by chromosome by group"""
    collapsed = collapsed[['chrom', 'density', 'extraction', 'digest', 'swga', 'rep', 'depth']]
    return collapsed.groupby(collapsed.columns.values[:-1].tolist()).\
        skew().\
        reset_index().\
        rename(columns={'depth' : 'skew'})
def wideform_chrom_pos(collapsed):
    """convert each chrom and position to a matrix of all samples"""
    wideform = pd.pivot_table(
        collapsed,
        index=['chrom', 'pos'],
        columns='sampleName',
        values='depth')
    return wideform
def plot_replicate_correlation(rep_corr):
    """plot replicate correlation of chromosomes as boxplot"""
    # keep only core chromosomes
    rep_corr = rep_corr[(rep_corr.chrom != 'Pf3D7_API_v3') & (rep_corr.chrom != 'Pf_M76611')]
    g = sns.catplot(
        x='density', y = 'correlation',
        kind='box', col='swga',
        row='extraction', height=10,
        hue='digest', data=rep_corr)
    g.savefig("../plots/replicate_correlation.png")
    plt.show()
    plt.close()
def plot_replicate_skewness(rep_skew):
    """plot skewness of chromosomes by group"""
    rep_skew = rep_skew[~((rep_skew.density == 100) & (rep_skew.extraction == 'Che') & (rep_skew.swga == 'Sof'))]
    rep_skew = rep_skew[(rep_skew.chrom != 'Pf3D7_API_v3') & (rep_skew.chrom != 'Pf_M76611')]

    # measure skewness
    g = sns.catplot(
        x='density', y = 'skew',
        col='swga', row='extraction',
        height=10, hue='digest',
        kind='box',
        data=rep_skew)
    g.savefig("../plots/sample_skewness.png")
    plt.show()
    plt.close()
def plot_chrom_coverage(wideform, core, chrom='Pf3D7_01_v3'):
    """plot chromosome depth as line plot across chromosome"""

    # normalize depth
    normed = (wideform / wideform.sum(axis=0) * 4e7).reset_index()
    normed = pd.melt(normed, id_vars=['chrom', 'pos'])

    # filter to chromosome
    core = core[core.chrom == chrom]
    normed = normed[normed.chrom == chrom]

    # plot median depth
    sns.lineplot(
        data=normed,
        x='pos',
        y='value',
        ci=95,
        n_boot=20,
        estimator=np.median)

    # add core intervals
    core.apply(lambda x : plt.hlines(-100, x.left_pos, x.right_pos), axis=1)

    # save figure
    plt.savefig("../plots/genome_cov/{0}.png".format(chrom))
    plt.show()


def main():
    collapsed_fn = "../data/opt4/mergedReplicates.collapsed.tab.gz"
    core_fn = "../resources/core.bed"
    core = pd.read_csv(core_fn, sep="\t", header=None, names=['chrom', 'left_pos', 'right_pos'])
    collapsed = pd.read_csv(collapsed_fn, sep='\t')
    collapsed = collapsed[(collapsed.chrom != 'Pf3D7_API_v3') & (collapsed.chrom != 'Pf_M76611')]

    collapsed_replicates = pivot_replicates(collapsed)
    rep_corr = replicate_correlation(collapsed_replicates)
    rep_skew = replicate_skew(collapsed)
    collapsed = add_sampleName(collapsed)

    wideform = wideform_chrom_pos(collapsed)

    plot_replicate_correlation(rep_corr)
    plot_replicate_skewness(rep_skew)
    [plot_chrom_coverage(wideform, core, chrom=i) for i in core.chrom.unique()]


if __name__ == '__main__':
    main()
