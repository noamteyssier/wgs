#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import vcf, sys
from tqdm import *
from scipy.spatial.distance import *

sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5}, style='ticks')


def iterVCF(fn):
    return vcf.Reader(open(fn, 'rb'))
def get_sampleAD(var):
    return np.array([s['AD'] for s in var.samples])
def get_sampleNames(fn):
    return np.array(vcf.Reader(open(fn, 'rb')).samples)
def get_singleSampleSnps(ad):
    snp_sample_idx = np.where(ad[:,1] > 0)[0]
    if snp_sample_idx.size == 1:
        return snp_sample_idx
    else:
        return None
def get_altFrequency(ad):
    return ad[:,1] / (ad[:,0] + ad[:,1])
def get_distance(maf_matrix):
    distances = pdist(maf_matrix.T, metric = lambda x,y : ((x>0)!=(y>0)).mean())
    return 1 - squareform(distances)
def get_af_frame(vcf_fn, sample_names):

    # single_sample_snp_count = np.zeros(sample_names.size)
    # single_sample_snp_idx = []
    alternative_frequencies = []

    # i = 0
    for v in tqdm(iterVCF(vcf_fn)):
        # sys.exit()

        sampleAD = get_sampleAD(v)
        # snp_sample_idx = get_singleSampleSnps(sampleAD)
        alt_frequencies = get_altFrequency(sampleAD)
        alternative_frequencies.append(alt_frequencies)

        # # count number of single sample snps
        # if snp_sample_idx:
        #     single_sample_snp_count[snp_sample_idx] += 1
        #     single_sample_snp_idx.append(i)
        # i+=1


    # single_sample_snp_idx = np.array(single_sample_snp_idx)
    # single_sample_snps = pd.DataFrame(single_sample_snp_count, index=sample_names)
    snp_frame = pd.DataFrame(alternative_frequencies, columns = sample_names)
    return snp_frame
def get_singleSampleSnps(snp_frame):
    print(snp_frame)
    sys.exit()

def main():
    vcf_fn = '../data/complete_set/optim.snps.recal.vcf.gz'
    ss_fn = '../data/complete_set/single_strains.snps.vcf.gz'
    error_rates = "../data/complete_set/snpConcordance.tab"
    errors = pd.read_csv(error_rates, sep="\t")


    np.unique(errors[['sample_i', 'sample_j']].values.reshape(-1))

    single_strains = errors.sample_i.unique()[-7:]
    samples = errors.sample_i.unique()[:-7]

    ss_distance = errors[errors.sample_i.isin(single_strains)]
    distance = ss_distance.pivot_table(index = 'sample_i', columns = 'sample_j', values='errorRate')


    sns.heatmap(1 - distance[distance.index.isin(['U51', 'fcr3'])].T)

    sample_distances = ss_distance[(ss_distance.sample_i.isin(['U51', 'fcr3'])) & (~ss_distance.sample_j.isin(single_strains))]
    max_similarity = 1 - sample_distances.groupby('sample_j').agg({'errorRate' : 'min'})
    sns.barplot(max_similarity.index, max_similarity.errorRate)


if __name__ == '__main__':
    main()
