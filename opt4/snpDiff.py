#!/usr/bin/env python3

import numpy as np
import pandas as pd
import vcf
import sys, os
# from scipy.spatial.distance import squareform, pdist
# from statsmodels.formula.api import ols
#
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

class SampleVariants:
    def __init__(self, dir, fn):
        self.dir = dir
        self.fn = fn
        self.sampleName = fn.strip('.vcf').strip('snpcore_')

        self.gt_conversion = {'A' : 0, 'C' : 1, 'T' : 2, 'G' : 3}

        self.open_vcf()
    def open_vcf(self):
        """read vcf into dataframe"""
        sample_df = []
        header = 0
        for record in vcf.Reader(open(self.dir + self.fn, 'r')):
            sample_df.append((self.sampleName, record.CHROM, record.POS, self.gt_conversion[str(record.ALT[0])], record.QUAL, record.INFO['AF1']))
            header += 1
            # if header == 1000:
            #     break
        self.vcf = pd.DataFrame(sample_df, columns = ['sampleName','chrom', 'pos', 'alt', 'qual', 'af'])
        self.vcf = self.vcf.set_index(['chrom', 'pos', 'sampleName'])
        # self.vcf.columns = col_index
def merge_variants(dir):
    vcf_fns = [i for i in os.listdir(dir) if '.vcf' in i]
    variant_calls = []
    for v in vcf_fns:
        print(v)
        sv = SampleVariants(dir, v)
        variant_calls.append(sv.vcf)

    merged_v = pd.concat(variant_calls).\
        reset_index()
    pivot_v = pd.pivot_table(merged_v, index = ['chrom', 'pos', 'sampleName'], values=['alt', 'af', 'qual'])
    pivot_v.to_csv("resources/merged_variants.tab.gz", sep="\t")
def plot_snp_sample_dist(merged_v):
    """distribution of number of samples across variants"""
    ns = merged_v.reset_index()[['chrom', 'pos', 'num_samples']].drop_duplicates()
    sns.distplot(ns.num_samples, kde=False)
    plt.show()
def single_sample_snps(merged_v):
    """visualization of single_sample_snps by group"""
    single_snps = merged_v.\
        groupby(level=2).\
        num_samples.\
        apply(lambda x : x[x==1].size).\
        reset_index().\
        rename(columns={'num_samples' : 'num_single_snps'})

    single_snps[['density', 'extraction', 'digest', 'swga', 'rep']] = \
        single_snps.apply(lambda x : x.sampleName.split('-'), axis = 1, result_type='expand')

    sns.catplot(
        data=single_snps,
        x = 'density',
        y = 'num_single_snps',
        col='swga',
        row='extraction',
        hue='digest',
        kind='bar',
        height=10)
    plt.show()
    plt.close()

def main():
    dir = "variant_calls/"
    # merge_variants(dir)
    merged_v = pd.read_csv("resources/merged_variants.tab.gz", sep='\t').set_index(['chrom', 'pos', 'sampleName'])

    # apply quality filter on snps
    merged_v = merged_v[merged_v.qual > 1]
    merged_v['num_samples'] = merged_v.groupby(level = [0,1]).apply(lambda x : x.shape[0])

    plot_snp_sample_dist(merged_v)
    single_sample_snps(merged_v)



if __name__ == '__main__':
    main()
