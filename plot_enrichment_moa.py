#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 16:24:11 2020

@author: jingxing
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import numpy as np

workdir = '/Users/jingxing/Documents/CoV/WARS/data/consensus_rank_200213/'
enrich_table = pd.read_csv(workdir + 'Enrichment_MoA_Targets.csv', index_col=0)
enrich_table['Num'] = [len(i.split('|')) for i in enrich_table.Median_Rank]
moa_table = enrich_table.loc[enrich_table.TYPE == 'MoA']
moa_table = moa_table.loc[moa_table.FDR < 0.25].iloc[:10]
moa_table = moa_table.sort_values(by='P', ascending=False)
tar_table = enrich_table.loc[(enrich_table.TYPE == 'Target') & (enrich_table.Num > 5)]
tar_table = tar_table.sort_values(by='P', ascending=False)


#plt.style.use('ggplot')
FIG = plt.figure(figsize=(5, 4), dpi=600)
x = -np.log10(moa_table.FDR.to_list())
#x = moa_table.FDR
y = np.arange(moa_table.shape[0])
size = moa_table.Num.to_list()
sct = plt.scatter(x, y, s=size)
leg = sct.legend_elements(prop="sizes", num=4, alpha=0.5)
plt.legend(*leg, title="Drug Number", loc='lower right', title_fontsize=6, fontsize=6)
plt.xlabel('-log10 FDR', fontsize=8)
plt.xticks(fontsize=6)
plt.yticks(ticks=y, labels=moa_table.index, fontsize=8)
plt.title('MoAs Enriched', fontsize=8)
FIG.tight_layout()
FIG.savefig(workdir + 'Plot_Enrichment_MoA.pdf', transparent=True)

'''
FIG, axs = plt.subplots(1, 2, figsize=(8, 6), dpi=300)
for i in range(2):
    axs[i].set_xlabel('-log10 P')
    axs[i].grid(True)
x = -np.log10(moa_table.P.to_list())
y = np.arange(moa_table.shape[0])
size = moa_table.Num.to_list()
sct = axs[0].scatter(x, y, s=size)
leg = sct.legend_elements(prop="sizes", num=4, alpha=0.5)
axs[0].legend(*leg, title="Drug Number", loc='lower right')
axs[0].set_yticks(y)
axs[0].set_yticklabels(moa_table.index)
axs[0].set_title('MoAs Enriched')
x = -np.log10(tar_table.P.to_list())
y = np.arange(tar_table.shape[0])
size = tar_table.Num.to_list()
sct = axs[1].scatter(x, y, s=size)
leg = sct.legend_elements(prop="sizes", num=4, alpha=0.5)
axs[1].legend(*leg, title="Drug Number", loc='lower right')
axs[1].set_yticks(y)
axs[1].set_yticklabels(tar_table.index)
axs[1].set_title('Targets Enriched')
FIG.tight_layout()
FIG.savefig(workdir + 'Plot_Enrichment_moa_target.pdf')
plt.show()
'''