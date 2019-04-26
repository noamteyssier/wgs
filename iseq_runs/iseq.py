#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})


class Balance:
    plates = []
    """rebalance an equivolume iseq plate for repooling"""
    def __init__(self, flagstat, min_mapping=0.1):
        self.flagstat = flagstat
        self.min_mapping = min_mapping

        self.__original_statistics__()
        self.__percentages__()
        self.__trim_undermapped__()
        self.__balance_mapped__()

        Balance.plates.append(self.flagstat)

    def __original_statistics__(self):
        """gather original statistics of the plate"""
        self.original_size = self.flagstat.shape[0]
        self.u_reads = self.flagstat[self.flagstat.sampleName=='Undetermined'].total.sum()
    def __filtered_statistics__(self):
        """gather statistics of the plate post processing"""
        pass
    def __percentages__(self):
        """percentage of mapped reads per sample"""
        self.flagstat['pc_mapped'] = self.flagstat.mapped / self.flagstat.total
    def __trim_undermapped__(self):
        """filter undermapped reads with threshold"""
        self.flagstat = self.flagstat[self.flagstat.pc_mapped >= self.min_mapping]
    def __balance_mapped__(self):
        """generate balancing coeff of all samples based off mapped read counts"""
        # ignore undetermined from total sum
        total_mapped = self.flagstat[self.flagstat.sampleName != 'Undetermined'].mapped.sum()
        self.flagstat['balance_coeff'] = 1 / (self.flagstat.mapped / total_mapped) * 0.05



def plate_summaries(flagstat):
    """flagstat summary by plate"""
    summ = flagstat.\
        groupby('plate').agg({
            'total' : ['sum', 'mean', 'std'],
            'mapped' : ['mean', 'std']})
    return summ

def main():
    flagstat_fn = "data/flagstat.tab"
    flagstat = pd.read_csv(flagstat_fn, sep="\t")

    summ = plate_summaries(flagstat)
    rebal = flagstat.\
        groupby('plate').\
        apply(lambda x : Balance(x))

    plates = pd.concat(Balance.plates)
    plates.to_csv(sys.stdout, sep='\t', index=False)

if __name__ == '__main__':
    main()
