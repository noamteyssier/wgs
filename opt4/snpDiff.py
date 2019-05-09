#!/usr/bin/env python3

import numpy as np
import pandas as pd
import vcf as pyvcf
import sys
from scipy.spatial.distance import squareform, pdist
from statsmodels.formula.api import ols

import seaborn as sns
import matplotlib.pyplot as plt
sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})

def vcf_mat(vcf):
    as_mat = []
    for read in vcf:
        arr = np.array([int(s['GT']) if s['GT'] != '.' else 0 for s in read.samples])
        as_mat.append(arr)

    mat = pd.DataFrame(np.vstack(as_mat), columns=vcf.samples).T
    return mat

def sample_dist(mat, plot=True):
    distances = pd.DataFrame(squareform(pdist(mat, 'sqeuclidean')))
    distances.index = distances.columns = vcf.samples
    if plot:
        sns.heatmap(distances)
        plt.show()
    return distances

def snp_sum(mat, plot=True):
    mat = mat.T
    mat['snp_sum'] = mat.apply(lambda x : x.sum(), axis = 1)
    if plot:
        sns.distplot(mat.snp_sum)
        plt.show()
    return mat

def sampleSingleSnps(mat):
    def single_snps(x):
        nonzero = np.nonzero(x.values)[0]
        if nonzero.size == 1:
            return nonzero[0]
    single_snp_samples = mat.\
        apply(lambda x : single_snps(x), axis = 0).\
        dropna()
    sample_vec = np.unique(single_snp_samples, return_counts = True)

    return pd.DataFrame(sample_vec[1].T, columns = ['num_singlesnps'], index=mat.index.values).\
        reset_index().\
        rename(columns={'index' : 'sampleName'})

def main():
    snp_fn = "data/snpDiff.tab"
    vcf_fn = "data/intersect_replicateVariance.vcf"
    vcf = pyvcf.Reader(open(vcf_fn, 'r'))

    # vcf.samples
    # vcf.infos

    mat = vcf_mat(vcf)
    distances = sample_dist(mat)
    ssum = snp_sum(mat)
    singleSnps = sampleSingleSnps(mat)

    singleSnps[['density','extraction','digest','swga','rep']] = singleSnps.apply(lambda x : x.sampleName.split('-'), result_type='expand', axis = 1)

    # single snp distribution by density
    [sns.distplot(singleSnps[singleSnps.density==i].num_singlesnps) for i in ['100', '1000']]

    # categorical differences with p1000
    sns.catplot(data=singleSnps[singleSnps.density=='1000'],
        row='extraction', col='swga', kind='violin', inner='stick',
        x='digest', y ='num_singlesnps', height=8)

    piv_singles = pd.pivot_table(singleSnps, index=['density', 'extraction', 'digest', 'swga'], columns='rep', values='num_singlesnps').\
        reset_index().\
        rename(columns={'1':'rep1', '2':'rep2'})

    # replicate differences
    piv_singles['rep_diff'] = np.abs(piv_singles['rep1'] - piv_singles['rep2'])
    [sns.distplot(piv_singles[(piv_singles.rep_diff <= 200) & (piv_singles.density == i)].rep_diff) for i in ['100', '1000']]
    [sns.distplot(piv_singles[(piv_singles.density == '1000') & (piv_singles.digest == i)].rep_diff, rug=True) for i in ['M', 'R']]


    # ANOVA testing of replicate differences
    anova_test = piv_singles[piv_singles.rep_diff < 300]
    ols('rep_diff ~ density + extraction + digest + swga', data= anova_test).fit().summary()

    # ANOVA testing of single_snps 
    ols('num_singlesnps ~ density + extraction + digest + swga', data=singleSnps).fit().summary()

if __name__ == '__main__':
    main()
