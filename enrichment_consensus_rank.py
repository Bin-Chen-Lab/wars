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

score_file = '/Users/jingxing/Documents/CoV/WARS/data/consensus_rank/consensus_rank_sRGES_drugs_CaseVSCtrl.csv'

drug_scores = pd.read_csv(score_file, index_col=0)
enr_type = 'moa'
score_type = 'Median_Rank'

drug_scores = drug_scores.sort_values(by=score_type)
drug_scores['RANK'] = range(1, drug_scores.shape[0] + 1)

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
    rank_a = drug_scores.loc[A, 'RANK'].astype(str)
    s, p = ranksums(score_a, score_b)
    p = 0.5 * p if s < 0 else 1 - 0.5 * p
    result.append([enr_set, p, '|'.join(rank_a)])
result = pd.DataFrame(data=result, columns=[enr_type, 'P', 'Ranks'])
result = result.sort_values(by='P')
result['FDR'] = multipletests([1 if np.isnan(i) else i for i in result.P], \
                              alpha=0.25, method='fdr_bh')[1]
result = result.loc[:, [enr_type, 'P', 'FDR', 'Ranks']]
result.to_csv('../data/consensus_rank/enrichment_%s.csv'%enr_type, index=False)
