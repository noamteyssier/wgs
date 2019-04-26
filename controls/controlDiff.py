#!/usr/bin/env python3

import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os, sys

from scipy.spatial.distance import pdist, squareform

sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})
pd.set_option('mode.chained_assignment', None) # remove settingwithcopywarning


def load_bcfs(dir):
    dir = dir.rstrip('/')
    fns = [f for f in os.listdir(dir) if 'SS' in f]
    return [pysam.VariantFile(dir + '/' + i) for i in fns], fns

def get_recs(bcf):
    vars = []
    for rec in bcf.fetch():
        if 'INDEL' in rec.info.keys():
            continue
        # if rec.chrom != 'Pf3D7_07_v3':
        #     continue

        depth = rec.info['DP']
        m_depth = np.array(rec.info['AC'])
        min_pc = (m_depth / depth)
        vars.append([rec.chrom, rec.pos, rec.ref, rec.alts[0], min_pc[0]])
    return pd.DataFrame(vars, columns = ['chrom', 'pos', 'ref', 'alts', 'min_pc'])


def main():
    dir = 'data/bcf'
    dir.rstrip('/')
    bcfs, fns = load_bcfs(dir)

    recs = []
    for i in range(len(bcfs)):
        rec = get_recs(bcfs[i])
        rec['sample'] = fns[i]
        recs.append(rec)

    # recs = [get_recs(b) for b in bcfs]

    df = pd.concat(recs)

    p = df.pivot_table(index=['chrom', 'pos'], columns = 'sample', values = 'min_pc').fillna(0)
    p = p.astype('bool').astype('int')
    dist = pd.DataFrame(squareform(pdist(p.T, 'jaccard')), columns = fns, index = fns)
    print(dist)

    sns.heatmap(1 - dist)
    # sns.heatmap(1 - dist / dist.max())
    plt.show()

if __name__ == '__main__':
    main()
