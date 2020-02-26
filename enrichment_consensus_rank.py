#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:30:20 2020

@author: jingxing
"""

import pandas as pd
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'


workdir = '../data/consensus_rank/'
score_file = workdir + 'consensus_rank_sRGES_drugs_CaseVSCtrl.csv'

drug_scores = pd.read_csv(score_file, index_col=0)
enr_type = 'moa'
score_type = 'Median_Rank'

drug_scores = drug_scores.sort_values(by=score_type)

tmp = drug_scores[enr_type].dropna().to_list()
enrich_sets = set()
for t in tmp:
    enrich_sets = enrich_sets.union(t.split('|'))
del tmp

result = []
for enr_set in enrich_sets:
    A = drug_scores[enr_type].str.contains(enr_set, na='')
    A = A.where(A == True).dropna().index
    score_a = drug_scores.loc[A, score_type]
    score_b = drug_scores.loc[~drug_scores.index.isin(A), score_type]
    rank_a = drug_scores.loc[A, score_type].astype(str)
    s, p = ranksums(score_a, score_b)
    p = 0.5 * p if s < 0 else 1 - 0.5 * p
    result.append([enr_set, p, '|'.join(rank_a)])
result = pd.DataFrame(data=result, columns=[enr_type, 'P', score_type])
result = result.sort_values(by='P')
result['FDR'] = multipletests([1 if np.isnan(i) else i for i in result.P], \
                              alpha=0.25, method='fdr_bh')[1]
result = result.loc[:, [enr_type, 'P', 'FDR', score_type]]
result.to_csv(workdir + 'enrichment_%s.csv'%enr_type, index=False)

# Keep top 10 enriched MoA
result['Num'] = [len(i.split('|')) for i in result[score_type]]
moa_table = result.loc[result.FDR < 0.25].iloc[:10]
moa_table = moa_table.sort_values(by='P', ascending=False)

# Plot
FIG = plt.figure(figsize=(5, 4), dpi=600)
x = -np.log10(moa_table.FDR.to_list())
y = np.arange(moa_table.shape[0])
size = moa_table.Num.to_list()
sct = plt.scatter(x, y, s=size)
leg = sct.legend_elements(prop="sizes", num=4, alpha=0.5)
plt.legend(*leg, title="Drug Number", loc='lower right', title_fontsize=6, fontsize=6)
plt.xlabel('-log10 FDR', fontsize=10)
plt.xticks(fontsize=6)
plt.yticks(ticks=y, labels=moa_table.index, fontsize=12)
plt.title('MoAs Enriched', fontsize=12)
FIG.tight_layout()
FIG.savefig(workdir + 'Fig2B_Enrichment_MoA.pdf', transparent=True)





