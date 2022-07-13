#!/usr/bin/env python
# coding: utf-8
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
## sgRNA_percentage
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

# 1. Check arguments
import argparse
parser = argparse.ArgumentParser(description="""
            The purpose of this script is to classify junction-spanning reads from STAR 2.7.9a
            alignment result. 
            SJ_out.tab file provide information of JSRs that support sgRNAs
            Basically, this script echo 5' and 3' breakpoint and reads count of the discontinous 
            transcription 
            """)
parser.add_argument("-c", "--canonical_count",required=True,
                    help="input canonical_sgRNA_count.tsv")
parser.add_argument("-s", "--noncanonical_count",required=True,
                    help="input nc_sgRNA_count.tsv")
parser.add_argument("-n", "--FileNamePrefix",required=True)

args = parser.parse_args()
canonical_count = args.canonical_count
noncanonical_count = args.noncanonical_count
FileNamePrefix = args.FileNamePrefix

#-------------------------------------------------------------------------------------------------#

# 2. Get input file

import pandas as pd
orfcolors = pd.read_csv('./colorcode-rainbow.txt', sep=' ', names=['orf', 'color'], index_col=0)['color'].to_dict()

canonical_counts = pd.read_csv(canonical_count, sep='\t', header=0 ,names=['orf', 'count'])

ORForder = list(orfcolors.keys())
ORForder = ['ORF1ab'] + ORForder[2:]
orfcolors['ORF1ab'] = orfcolors['ORF1a']
canonical_counts['order'] = canonical_counts['orf'].apply(lambda x: ORForder.index(x))
canonical_counts['color'] = canonical_counts['orf'].apply(lambda x: orfcolors[x])

nonproductive_counts = pd.read_csv(noncanonical_count, sep='\t')
import numpy as np
nonproductive_counts['order'] = np.arange(len(nonproductive_counts)) + 100
nonproductive_counts['color'] = ['#' + ('%02x' % int(i) * 3)
                                 for i in np.linspace(0, 200, len(nonproductive_counts))]

#-------------------------------------------------------------------------------------------------#

# 3. Combine canonical and noncanonical sgRNA into one dataframe

fullcounts = pd.concat([canonical_counts, nonproductive_counts]).reset_index(drop=True)
fullcounts['count'].sum()

fullcounts['pct'] = fullcounts['count'] / fullcounts['count'].sum() * 100

fullcounts['bottom'] = [0] + fullcounts['pct'].cumsum().tolist()[:-1]
fullcounts['top'] = fullcounts['bottom'] + fullcounts['pct']

fullcounts = fullcounts.sort_values(by='order').copy()

#-------------------------------------------------------------------------------------------------#

# 4. Draw final plot

from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt

fig, ax = plt.subplots(1, 1, figsize=(5.6, 1.5))

for _, row in fullcounts.iterrows():
    ax.barh([0.5], [row['pct']], 1, left=[row['bottom']],
           fc=row['color'], zorder=1)
    if  row['pct'] >=10:
        ax.annotate(row['orf'].replace('ORF', ''), ((row['bottom'] + row['top'])/2, 0.5),
                    ha='center', fontsize=12, va='center')
#row['order'] <=9 
ax.set_xlim(0, 100)
for spname in 'top left right bottom'.split():
    ax.spines[spname].set_visible(False)

plt.setp(ax.get_yticklines(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_xlabel('Junction-spanning reads (%)')

plt.tight_layout()
plt.savefig(FileNamePrefix + "sgRNA_precentage.jpg")


