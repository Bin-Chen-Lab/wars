#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:52:21 2020

@author: jingxing
"""

import pandas as pd
import os
from scipy.stats import ranksums, spearmanr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns


pos_ctrl_ec50 = pd.read_csv('../data/raw/pos_ctrl_EC50.tsv', sep='\t', index_col=0)
pos_ctrl_ec50['log10_EC50'] = np.log10(pos_ctrl_ec50.EC50 * 10**(-6))
pos_ctrls = pos_ctrl_ec50.index.to_list()
inpdir = '../data/all_comparisons/'
drug_info_file = '../data/raw/repurposing_drugs_20180907.txt'
outdir = '../data/consensus_rank_200213/'
if not os.path.exists(outdir): os.mkdir(outdir)
encd = 'mac_roman'

# Get the sRGES ranks and scores of the positive control drugs for all the infection comparison signatures
num_dz_g, pct_top, drug_lst = 50, 0.05, None
rank_dict, score_dict = dict(), dict()
for comp in os.listdir(inpdir):
    if all([not comp.startswith(st) for st in ['all_virus_', 'GSE', 'meta_GSE']]): continue
    drug_score_file = inpdir + comp + '/sRGES_drugs.csv'
    if not os.path.exists(drug_score_file): continue
    score_lst = pd.read_csv(drug_score_file, encoding=encd)
    score_lst = score_lst.loc[score_lst.n > 1]
    #score_lst.name = score_lst.name.str.lower()
    if not drug_lst: drug_lst = sorted(score_lst.name.to_list())
    score_lst.index = range(1, score_lst.shape[0] + 1)
    drug_scores = score_lst.loc[score_lst.name.isin(pos_ctrls), ['name', 'sRGES']]
    rank_dict[comp] = dict(zip(drug_scores.name, drug_scores.index))
    score_dict[comp] = dict(zip(drug_scores.name, drug_scores.sRGES))
pos_ctrl_ranks = []
pos_ctrls = [d for d in pos_ctrls if d in drug_lst]
pos_ctrl_ec50 = pos_ctrl_ec50.loc[pos_ctrls]
for k in rank_dict:
    dzsig_file = inpdir+k+'/dz_signature_lincs_used.csv'
    if not os.path.exists(dzsig_file): continue
    dzsig = pd.read_csv(dzsig_file, encoding=encd)
    pos_ctrl_ranks.append([k, dzsig.shape[0]] + [rank_dict[k][d] for d in pos_ctrls])
pos_ctrl_ranks = pd.DataFrame(data=pos_ctrl_ranks, columns=['model', 'Num_Genes'] + pos_ctrls)
# sRGES of the positive drugs across all the comparisons
pos_ctrl_scores = []
pos_ctrls_1 = pos_ctrl_ec50.dropna().index.to_list()
for k in score_dict:
    pos_ctrl_scores.append([k] + [score_dict[k][d] for d in pos_ctrls_1])
pos_ctrl_scores = pd.DataFrame(data=pos_ctrl_scores, columns=['model'] + pos_ctrls_1)
pos_ctrl_scores.index = pos_ctrl_scores.model
pos_ctrl_scores = pos_ctrl_scores.iloc[:, 1:]
pos_ctrl_scores.to_csv(outdir + 'sRGES_positiveCtrls.csv')

# Select the comparisons with enough dysregulated genes
pos_ctrl_ranks = pos_ctrl_ranks.loc[pos_ctrl_ranks.Num_Genes >= num_dz_g]
pos_ctrl_ranks.index = pos_ctrl_ranks.model
pos_ctrl_scores = pos_ctrl_scores.loc[pos_ctrl_scores.index.isin(pos_ctrl_ranks.index)]

# Select the comparisons can enrich positive control drugs
tmp = pos_ctrl_ranks.loc[:, pos_ctrls]
enrich_test = []
for comp in pos_ctrl_ranks.index:
    pos = tmp.loc[comp].to_list()
    back = [r for r in range(1, len(drug_lst) + 1) if r not in pos]
    s, p = ranksums(pos, back)
    p = p * 0.5 if s < 0 else 1 - p * 0.5
    enrich_test.append(p)
del tmp
pos_ctrl_ranks['P_Enrichment'] = enrich_test
# Calculate Spearman r of srges vs ec50 within each comparison
cor_r, cor_p = [], []
for comp in pos_ctrl_scores.index:
    scores = pos_ctrl_scores.loc[comp].to_list()
    ec50 = pos_ctrl_ec50.loc[pos_ctrl_scores.columns, 'log10_EC50']
    sr, p = spearmanr(scores, ec50)
    p = p * 0.5 if sr > 0 else 1 - p * 0.5
    cor_r.append(sr)
    cor_p.append(p)
pos_ctrl_ranks['SpearmanR_sRGES_EC50'] = cor_r
pos_ctrl_ranks['P_SpearmanR'] = cor_p
pos_ctrl_ranks = pos_ctrl_ranks.sort_values(by='SpearmanR_sRGES_EC50', ascending=False)
pos_ctrl_ranks.to_csv(outdir+'comparisons_validation.csv', index=False)

cmpr_sele = pos_ctrl_ranks.loc[(pos_ctrl_ranks['P_Enrichment'] < 0.05) & \
                               (pos_ctrl_ranks['P_SpearmanR'] < 0.05) & \
                               (pos_ctrl_ranks['SpearmanR_sRGES_EC50'] > 0.4)].index.to_list()
pos_ctrl_ranks.loc[cmpr_sele].to_csv(outdir+'comparisons_selected.csv')

# Plot the enrichment density and correlation score vs EC50
FIG = plt.figure(figsize=(10, len(cmpr_sele) * 4), dpi=300)
gs = gridspec.GridSpec(len(cmpr_sele) * 2, 2, \
                       height_ratios=[2 if i%2==0 else 1 for i in range(len(cmpr_sele) * 2)])
for i, comp in enumerate(cmpr_sele):
    # Draw density of positive drugs
    ax_dist = FIG.add_subplot(gs[i*2, 0])
    p = pos_ctrl_ranks.loc[comp, 'P_Enrichment']
    sns.distplot(pos_ctrl_ranks.loc[comp, pos_ctrls], hist=False, \
                 kde_kws={'kernel': 'triw'}, ax=ax_dist, \
                 label='P = %.2E'%p)
    ax_dist.set_xlim(0, len(drug_lst))
    ax_dist.set_axis_off()
    ax_dist.set_title(comp)
    # Draw rank positions of positive drugs
    ax_bcd = FIG.add_subplot(gs[i*2 + 1, 0])
    ax_bcd.eventplot(pos_ctrl_ranks.loc[comp, pos_ctrls])
    ax_bcd.set_xlim(0, len(drug_lst))
    ax_bcd.set_ylim(0.5, 1)
    ax_bcd.set_xlabel('Positive Drugs Rank')
    ax_bcd.get_yaxis().set_visible(False)
    for spine in ["left", "top", "right"]: ax_bcd.spines[spine].set_visible(False)
    # Draw correlation of srges and EC50
    ax_cor = FIG.add_subplot(gs[i*2:(i+1)*2, 1])
    x = pos_ctrl_scores.loc[comp]
    y = pos_ctrl_ec50.loc[pos_ctrl_scores.columns, 'log10_EC50']
    corr = pos_ctrl_ranks.loc[comp, 'SpearmanR_sRGES_EC50']
    p = pos_ctrl_ranks.loc[comp, 'P_SpearmanR']
    sns.regplot(x, y, label='Spearman R = %.2f\nP = %.4f'%(corr, p), ax=ax_cor)
    ax_cor.set_xlabel('sRGES')
    ax_cor.set_ylabel('log10 EC50')
    ax_cor.set_title(comp)
    ax_cor.legend()
FIG.tight_layout()
FIG.savefig(outdir + 'Plot_Enrichment_EC50_sRGES.pdf')
plt.close()


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


# Generate consensus rank table of drug sRGES on negative control signatures (MOCK)
rank_dict, cmpr_neg_ctrl = dict(), []
for comp in os.listdir(inpdir):
    if not comp.startswith('all_mock'): continue
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
conss_rank_n['Median_Rank_MOCK'] = conss_rank_n.median(axis=1)
conss_rank_n.index = [name.lower() for name in conss_rank_n.name]
conss_rank_n = pd.concat([conss_rank_n, drug_info], axis=1).loc[conss_rank_n.index]
conss_rank_n.to_csv(outdir + 'consensus_rank_sRGES_drugs_MOCK.csv', index=False)


# Generate consensus rank table of drug sRGES on other control signatures H1N1
rank_dict, cmpr_ctrl1, DZ = dict(), [], 'H1N1'
for comp in os.listdir(inpdir):
    if DZ not in comp.upper(): continue
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
    cmpr_ctrl1.append(comp)
conss_rank_n1 = []
for drug in drug_lst:
    conss_rank_n1.append([drug] + [rank_dict[comp][drug] for comp in cmpr_ctrl1])
conss_rank_n1 = pd.DataFrame(data=conss_rank_n1, columns=['name'] + cmpr_ctrl1)
conss_rank_n1['Median_Rank_'+DZ] = conss_rank_n1.median(axis=1)
conss_rank_n1.index = [name.lower() for name in conss_rank_n1.name]
conss_rank_n1 = pd.concat([conss_rank_n1, drug_info], axis=1).loc[conss_rank_n1.index]
conss_rank_n1.to_csv(outdir + 'consensus_rank_sRGES_drugs_%s.csv'%DZ, index=False)

# Combine activity table and toxicity/selectivity table
col_a = ['name', 'Median_Rank']
col_b = ['Median_Rank_MOCK']
col_c = ['Median_Rank_'+DZ] + drug_info.columns.to_list()
com = pd.concat([conss_rank.loc[:, col_a], conss_rank_n.loc[:, col_b], conss_rank_n1.loc[:, col_c]], axis=1)
com.to_csv(outdir + 'consensus_rank_sRGES_drugs_CaseVSCtrl.csv', index=False)









