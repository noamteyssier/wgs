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
    sns.catplot(
        x='density', y = 'correlation',
        kind='box', col='swga',
        row='extraction', height=10,
        hue='digest', data=rep_corr)
    plt.show()
    plt.close()
def plot_replicate_skewness(rep_skew):
    """plot skewness of chromosomes by group"""
    rep_skew = rep_skew[~((rep_skew.density == 100) & (rep_skew.extraction == 'Che') & (rep_skew.swga == 'Sof'))]
    rep_skew = rep_skew[(rep_skew.chrom != 'Pf3D7_API_v3') & (rep_skew.chrom != 'Pf_M76611')]

    # measure skewness
    sns.catplot(
        x='density', y = 'skew',
        col='swga', row='extraction',
        height=10, hue='digest',
        kind='box',
        data=rep_skew)
    plt.show()
    plt.close()
def plot_chrom_coverage(wideform):
    """plot chromosome coverage as line plot across genome"""
    sample_means = wideform.apply(lambda x : x.mean(), axis=0)
    normed = wideform * sample_means
    normed[['pos_mean', 'pos_var']] = normed.apply(lambda x : [x.mean(), x.var()], axis=1, result_type='expand')


    longform = normed[['pos_mean', 'pos_var']].reset_index().melt(id_vars=['chrom', 'pos'], var_name='summary_stat')
    longform['value'] = np.log10(longform['value'] + 1)
    g = sns.FacetGrid(longform, row='chrom', col='summary_stat',
        aspect=5, height=3, sharex=False, sharey=False,
        hue='summary_stat')
    g.map(sns.lineplot, "pos", "value")
    plt.show()
    plt.close()



def main():
    collapsed_fn = "data/mergedReplicates.collapsed.tab.gz"
    collapsed = pd.read_csv(collapsed_fn, sep='\t')
    collapsed = collapsed[(collapsed.chrom != 'Pf3D7_API_v3') & (collapsed.chrom != 'Pf_M76611')]

    collapsed_replicates = pivot_replicates(collapsed)
    rep_corr = replicate_correlation(collapsed_replicates)
    rep_skew = replicate_skew(collapsed)
    collapsed = add_sampleName(collapsed)

    wideform = wideform_chrom_pos(collapsed)
    wideform_100 = wideform.loc[:,[i for i in wideform.columns.values if '100-' in i]]
    wideform_1000 = wideform.loc[:,[i for i in wideform.columns.values if '1000-' in i]]

    plot_replicate_correlation(rep_corr)
    plot_replicate_skewness(rep_skew)
    plot_chrom_coverage(wideform)
    plot_chrom_coverage(wideform_100)
    plot_chrom_coverage(wideform_1000)




if __name__ == '__main__':
    main()
