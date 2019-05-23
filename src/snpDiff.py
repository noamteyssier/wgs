#!/usr/bin/env python3

import numpy as np
import pandas as pd
import vcf, sys, os
from scipy.spatial.distance import *

import seaborn as sns
import matplotlib.pyplot as plt
sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5}, style='ticks')

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
def merge_variants(dir, ofn="../resources/opt4_merged_variants.tab.gz"):
    """iterate through directory and concatenate vcfs to a single dataframe"""
    vcf_fns = [i for i in os.listdir(dir) if '.vcf' in i]
    variant_calls = []
    for v in vcf_fns:
        print(v)
        sv = SampleVariants(dir, v)
        variant_calls.append(sv.vcf)

    merged_v = pd.concat(variant_calls).\
        reset_index()
    pivot_v = pd.pivot_table(merged_v, index = ['chrom', 'pos', 'sampleName'], values=['alt', 'af', 'qual'])
    pivot_v.to_csv(ofn, sep="\t")
def plot_snp_sample_dist(merged_v):
    """distribution of number of samples across variants"""
    ns = merged_v.reset_index()[['chrom', 'pos', 'num_samples']].drop_duplicates()
    g = sns.distplot(ns.num_samples, kde=False)
    g.get_figure().savefig("../plots/snp_sample_numbers.png", height=8, width=8)
    plt.show()
def single_sample_snps(merged_v):
    """visualization of single_sample_snps by group"""
    single_snps = merged_v.\
        groupby('sampleName').\
        num_samples.\
        apply(lambda x : x[x==1].size).\
        reset_index().\
        rename(columns={'num_samples' : 'num_single_snps'})

    single_snps[['density', 'extraction', 'digest', 'swga', 'rep']] = \
        single_snps.apply(lambda x : x.sampleName.split('-'), axis = 1, result_type='expand')

    g = sns.catplot(
        data=single_snps,
        x = 'density',
        y = 'num_single_snps',
        col='swga',
        row='extraction',
        hue='digest',
        kind='bar',
        height=10)
    g.savefig("../plots/single_sample_snps.png")
    plt.show()
    plt.close()
def replicate_concordance(merged_v):
    """SNP concordance (measured as 1 - Jaccard distance) between replicates"""
    def jacc_sim(x):
        return(1 - jaccard(x.af['1'], x.af['2']))
    mv = merged_v.copy()
    mv[['density', 'extraction', 'digest', 'swga', 'rep']] = mv.apply(lambda x : x.sampleName.split('-'), axis = 1, result_type='expand')
    rep_mv = pd.pivot_table(mv, index=['chrom', 'pos', 'density', 'extraction', 'digest', 'swga'], columns='rep', values=['af']).reset_index().fillna(0)

    snp_concordance = rep_mv.\
        groupby(['chrom', 'density', 'extraction', 'digest', 'swga']).\
        apply(lambda x : jacc_sim(x)).\
        reset_index().\
        rename(columns = {0 : 'concordance'})

    g = sns.catplot(
        data=snp_concordance,
        col='density',
        row='extraction',
        hue='digest',
        x='swga',
        y='concordance',
        kind='boxen',
        height=10)
    g.set_titles("{row_name} + {col_name}p")
    g.set_ylabels("1 - Jaccard Distance")
    g.set_xlabels("sWGA Method")
    g.savefig("../plots/replicate_snp_concordance.png")
    plt.show()
    plt.close()
def single_strain_overlap(sample, unique_ss):
    """calculate percentage overlap with each single strain"""
    exp = sample[['chrom', 'pos', 'af']]
    return unique_ss.groupby('sampleName').\
        apply(lambda x : exp.merge(
            x[['chrom', 'pos', 'af']],
            how='inner',
            on=['chrom', 'pos']
            ).shape[0] / x.shape[0])
def quality_filter(snp_df, threshold=100):
    """basic quality filter for snps"""
    snp_df = snp_df[snp_df.qual > threshold]
    return snp_df
def count_samples(snp_df):
    """group by chrom & pos to count number of samples"""
    num_samples = snp_df.\
        groupby(['chrom', 'pos']).\
        sampleName.\
        count().\
        reset_index().\
        rename(columns = {'sampleName' : 'num_samples'})
    snp_df = snp_df.merge(num_samples, how = 'left')
    return snp_df
def preprocess_variant_df(snp_df, known_snps, threshold=100):
    """apply quality filtering and sample count"""
    snp_df = quality_filter(snp_df, threshold = threshold)
    snp_df = snp_df.\
        merge(known_snps[['chrom', 'pos', 'qual']], on=['chrom', 'pos'], how='left').\
        dropna()
    snp_df = count_samples(snp_df)
    return snp_df
def preprocess_known_snps(known_snps, threshold=1e4):
    """apply quality treshold to known snps and return"""
    known_snps = known_snps[known_snps.qual >= threshold]
    return known_snps
def snp_concordance(sdf1, sdf2):
    """calculate single strain concordance for all single strains for all samples"""
    ssc = sdf1.\
        groupby('sampleName').\
        apply(lambda x : single_strain_overlap(x, sdf2)).\
        reset_index()
    del ssc.columns.name

    return ssc
def optim_concordance_heatmap(v_ssc):
    """draw heatmap of snp concordance across optim"""
    df = v_ssc.set_index('sampleName')
    sns.heatmap(
        data=df,
        robust=True,
        square=True,
        center=0.7,
        linewidths=0.25,
        xticklabels=False)
    plt.ylabel("Sample Name")
    plt.savefig("../plots/grouped_snp_concordnace.png")



def main():
    exp_dir = "../variant_calls/opt4/"
    ss_dir = "../variant_calls/single_strains/"
    cc_dir = "../variant_calls/combinations/"

    # merge_variants(exp_dir)
    # merge_variants(ss_dir, ofn="../resources/ss_merged_variants.tab.gz")
    # merge_variants(cc_dir, ofn="../resources/cc_merged_variants.tab.gz")

    merged_v = pd.read_csv("../resources/opt4_merged_variants.tab.gz", sep='\t')
    merged_ss = pd.read_csv("../resources/ss_merged_variants.tab.gz", sep="\t")
    merged_cc = pd.read_csv("../resources/cc_merged_variants.tab.gz", sep="\t")
    known_snps = pd.read_csv("../resources/pf3k_snps.tab.gz", sep="\t")

    known_snps = preprocess_known_snps(known_snps)
    merged_v = preprocess_variant_df(merged_v, known_snps, threshold=200)
    merged_ss = preprocess_variant_df(merged_ss, known_snps, threshold=200)
    merged_cc = preprocess_variant_df(merged_cc, known_snps, threshold=200)

    unique_ss = merged_ss[merged_ss.num_samples == 1]
    ssc = snp_concordance(merged_v, unique_ss)
    v_ssc = snp_concordance(merged_v, merged_v)
    ss_ssc = snp_concordance(merged_ss, merged_ss)
    cc_ssc = snp_concordance(merged_cc, unique_ss)

    # distribution of overlaps with single strains
    [sns.distplot(np.log10(1e-6 + ssc.iloc[:,i])) for i in range(1, ssc.shape[1])]

    # replicate + sample concordance
    replicate_concordance(merged_v)
    optim_concordance_heatmap(v_ssc)

    # # plot_snp_sample_dist(merged_v)
    # # single_sample_snps(merged_v)


if __name__ == '__main__':
    main()
