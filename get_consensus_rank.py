#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:52:21 2020

@author: jingxing
"""

import pandas as pd
import os

pos_ctrls = ['niclosamide', 'lopinavir', 'ritonavir']
inpdir = '../data/comparisons/'
drug_info_file = '/Users/jingxing/Documents/CoV/WARS/data/raw/repurposing_drugs_20180907.txt'
outdir = '../data/consensus_rank/'
if not os.path.exists(outdir): os.mkdir(outdir)
encd = 'mac_roman'

# Get the sRGES ranks of the positive control drugs for all the infection comparison signatures
num_dz_g, pct_top, drug_lst = 50, 0.05, None
rank_dict = dict()
for comp in os.listdir(inpdir):
    if 'MOCK' in comp: continue
    drug_score_file = inpdir + comp + '/sRGES_drugs.csv'
    if not os.path.exists(drug_score_file): continue
    score_lst = pd.read_csv(drug_score_file, encoding=encd)
    score_lst = score_lst.loc[score_lst.n > 1]
    if not drug_lst: drug_lst = sorted(score_lst.name.to_list())
    score_lst.index = range(1, score_lst.shape[0] + 1)
    drug_scores = score_lst.loc[score_lst.name.isin(pos_ctrls), ['name', 'sRGES']]
    rank_dict[comp] = dict(zip(drug_scores.name, drug_scores.index))
pos_ctrl_ranks = []
for k in rank_dict:
    dzsig_file = inpdir+k+'/dz_signature_lincs_used.csv'
    if not os.path.exists(dzsig_file): continue
    dzsig = pd.read_csv(dzsig_file, encoding=encd)
    pos_ctrl_ranks.append([k, dzsig.shape[0]] + [rank_dict[k][d] for d in pos_ctrls])
pos_ctrl_ranks = pd.DataFrame(data=pos_ctrl_ranks, columns=['model', 'Num_Genes'] + pos_ctrls)

# Select the comparisons with enough dysregulated genes
pos_ctrl_ranks = pos_ctrl_ranks.loc[pos_ctrl_ranks.Num_Genes >= num_dz_g]
pos_ctrl_ranks.index = pos_ctrl_ranks.model
# Select the comparisons can rank any positive control drugs top 5%
tmp = pos_ctrl_ranks.loc[:, pos_ctrls]
cmpr_sele = tmp.where(tmp < len(drug_lst) * pct_top).dropna(how='all').index.to_list()
output = pos_ctrl_ranks.loc[cmpr_sele]
output.to_csv(outdir+'comparisons_selected.csv')
del tmp

# Generate consensus rank table containing all drugs
rank_dict = dict()
for comp in cmpr_sele:
    drug_score_file = inpdir + comp + '/sRGES_drugs.csv'
    score_lst = pd.read_csv(drug_score_file, encoding=encd)
    score_lst = score_lst.loc[score_lst.n > 1]
    score_lst.index = range(1, score_lst.shape[0] + 1)
    rank_dict[comp] = dict(zip(score_lst.name, score_lst.index))
conss_rank = []
for drug in drug_lst:
    conss_rank.append([drug] + [rank_dict[comp][drug] for comp in cmpr_sele])
conss_rank = pd.DataFrame(data=conss_rank, columns=['name'] + cmpr_sele)
conss_rank['Median_Rank'] = conss_rank.median(axis=1)
conss_rank['Q25_Rank'] = conss_rank.loc[:, cmpr_sele].quantile(q=0.25, axis=1)
conss_rank['Q75_Rank'] = conss_rank.loc[:, cmpr_sele].quantile(q=0.75, axis=1)
conss_rank.index = [name.lower() for name in conss_rank.name]

# Map drug information to the consensus rank table
drug_info = [l.strip().split('\t') for l in open(drug_info_file, encoding=encd) if not l.startswith('!')]
drug_info = pd.DataFrame(data=drug_info[1:], columns=drug_info[0])
drug_info.index = [name.lower().replace('"','') for name in drug_info.pert_iname]
drug_info = drug_info.iloc[:, 1:]
for c in drug_info.columns:
    drug_info[c] = drug_info[c].str.replace('"', '')
conss_rank = pd.concat([conss_rank, drug_info], axis=1).loc[conss_rank.index]
conss_rank.to_csv(outdir + 'consensus_rank_sRGES_drugs.csv', index=False)


# Generate consensus rank table of drug sRGES on negative control signatures
rank_dict, cmpr_neg_ctrl = dict(), []
for comp in os.listdir(inpdir):
    if 'MOCK' not in comp: continue
    dzsig_file = inpdir + comp +'/dz_signature_lincs_used.csv'
    if not os.path.exists(dzsig_file): continue
    dzsig = pd.read_csv(dzsig_file, encoding=encd)
    if dzsig.shape[0] < num_dz_g: continue
    drug_score_file = inpdir + comp + '/sRGES_drugs.csv'
    if not os.path.exists(drug_score_file): continue
    score_lst = pd.read_csv(drug_score_file, encoding=encd)
    score_lst = score_lst.loc[score_lst.n > 1]
    score_lst.index = range(1, score_lst.shape[0] + 1)
    rank_dict[comp] = dict(zip(score_lst.name, score_lst.index))
    cmpr_neg_ctrl.append(comp)
conss_rank_n = []
for drug in drug_lst:
    conss_rank_n.append([drug] + [rank_dict[comp][drug] for comp in cmpr_neg_ctrl])
conss_rank_n = pd.DataFrame(data=conss_rank_n, columns=['name'] + cmpr_neg_ctrl)
conss_rank_n['Median_Rank_NegCtrl'] = conss_rank_n.median(axis=1)
conss_rank_n['Q25_Rank_NegCtrl'] = conss_rank_n.loc[:, cmpr_neg_ctrl].quantile(q=0.25, axis=1)
conss_rank_n['Q75_Rank_NegCtrl'] = conss_rank_n.loc[:, cmpr_neg_ctrl].quantile(q=0.75, axis=1)
conss_rank_n.index = [name.lower() for name in conss_rank_n.name]
conss_rank_n = pd.concat([conss_rank_n, drug_info], axis=1).loc[conss_rank_n.index]
conss_rank_n.to_csv(outdir + 'consensus_rank_sRGES_drugs_negCtrl.csv', index=False)

# Combine activity table and toxicity table
col_a = ['name', 'Median_Rank', 'Q25_Rank', 'Q75_Rank']
col_b = ['Median_Rank_NegCtrl', 'Q25_Rank_NegCtrl', 'Q75_Rank_NegCtrl'] + drug_info.columns.to_list()
com = pd.concat([conss_rank.loc[:, col_a], conss_rank_n.loc[:, col_b]], axis=1)
com.to_csv(outdir + 'consensus_rank_sRGES_drugs_CaseVSCtrl.csv', index=False)









