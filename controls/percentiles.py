#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os, sys
sns.set(rc={'figure.figsize':(9, 6), 'lines.linewidth': 5})


def arrange_metrics(metrics):
    """convert wideform metrtics to longform"""
    cols = [c for c in metrics.columns if 'PCT_' in c and c[-1] == 'X'] + ['sample']
    p = pd.melt(metrics.loc[:,cols], id_vars='sample', var_name='depth', value_name='pc_genome')
    p['depth'] = p.depth.apply(lambda x : x.replace('PCT_', '').replace('X', '')).astype('int')
    return p
def split_samplename(metrics):
    """split samplename into combination and density"""
    def sample_split(x):
        return x['sample'].split('-')
    metrics[['combination', 'density']]= metrics.apply(
        sample_split, axis = 1, result_type='expand')
    return metrics

def main():
    metrics = pd.read_csv("data/percentilesCoverage.tab", sep="\t")
    long_form = arrange_metrics(metrics)
    long_form = split_samplename(long_form)

    sns.lineplot(data=long_form, x='depth', y='pc_genome', style='density')
    plt.xticks(np.linspace(0, 100, 21))
    plt.savefig("plots/percentiles.pdf")

    long_form.to_csv("data/control_percentiles.tab", sep="\t", index=False)

if __name__ == '__main__':
    main()
