#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})

def main():
    df = pd.read_csv("data/moz_p1_FA.tab", sep="\t")

    sns.lineplot(df, x= 'ng.uL', y ='')
    pass

if __name__ == '__main__':
    main()
