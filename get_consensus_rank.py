#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:52:21 2020

@author: jingxing
"""

import pandas as pd
import os
from scipy.stats import ranksums, spearmanr
from statsmodels.stats.multitest import multipletests
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
import seaborn as sns
import re


pos_ctrl_ec50 = pd.read_csv('../data/raw/pos_ctrl_EC50.tsv', sep='\t', index_col=0)
pos_ctrl_ec50['log10_EC50'] = np.log10(pos_ctrl_ec50.EC50 * 10**(-6))
pos_ctrls = pos_ctrl_ec50.index.to_list()
inpdir = '../data/all_comparisons/'
drug_info_file = '../data/raw/repurposing_drugs_20180907.txt'
outdir = '../data/consensus_rank/'
if not os.path.exists(outdir): os.mkdir(outdir)
encd = 'mac_roman'

# Get the sRGES ranks and scores of the positive control drugs for all the infection comparison signatures
num_dz_g, num_comp, drug_lst = 50, 0, None
rank_dict, score_dict = dict(), dict()
for comp in os.listdir(inpdir):
    drug_score_file = inpdir + comp + '/sRGES_drugs.csv'
    if not os.path.exists(drug_score_file): continue
    score_lst = pd.read_csv(drug_score_file, encoding=encd)
    score_lst = score_lst.loc[score_lst.n > 1]
    if score_lst.sRGES.iloc[0] > -0.25: continue
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
pos_ctrl_ranks.index = pos_ctrl_ranks.model
pos_ctrl_ranks = pos_ctrl_ranks.iloc[:, 1:]
# sRGES of the positive drugs across all the comparisons
pos_ctrl_scores = []
pos_ctrls_1 = pos_ctrl_ec50.dropna().index.to_list()
for k in score_dict:
    pos_ctrl_scores.append([k] + [score_dict[k][d] for d in pos_ctrls_1])
pos_ctrl_scores = pd.DataFrame(data=pos_ctrl_scores, columns=['model'] + pos_ctrls_1)
pos_ctrl_scores.index = pos_ctrl_scores.model
pos_ctrl_scores = pos_ctrl_scores.iloc[:, 1:]
pos_ctrl_scores.to_csv(outdir + 'sRGES_positiveCtrls.csv')

# Select the comparisons can enrich positive control drugs
tmp = pos_ctrl_ranks.loc[:, pos_ctrls]
enrich_test = []
for comp in pos_ctrl_ranks.index:
    # Must include enough DE genes
    if pos_ctrl_ranks.loc[comp, 'Num_Genes'] < num_dz_g:
        enrich_test.append(1.0)
        continue
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
    # Must include enough DE genes
    if pos_ctrl_ranks.loc[comp, 'Num_Genes'] < num_dz_g:
        cor_r.append(0.0)
        cor_p.append(1.0)
        continue
    scores = pos_ctrl_scores.loc[comp].to_list()
    ec50 = pos_ctrl_ec50.loc[pos_ctrl_scores.columns, 'log10_EC50']
    sr, p = spearmanr(scores, ec50)
    p = p * 0.5 if sr > 0 else 1 - p * 0.5
    cor_r.append(sr)
    cor_p.append(p)
pos_ctrl_ranks['SpearmanR_sRGES_EC50'] = cor_r
pos_ctrl_ranks['P_SpearmanR'] = cor_p
pos_ctrl_ranks = pos_ctrl_ranks.sort_values(by='SpearmanR_sRGES_EC50', ascending=False)

cov_posctrl_ranks = pos_ctrl_ranks.loc[((pos_ctrl_ranks.index.str.contains('SARS')) | \
                                       (pos_ctrl_ranks.index.str.contains('MERS'))) & \
                                       (~pos_ctrl_ranks.index.str.startswith('all_mock_'))]
cov_posctrl_ranks['FDR_Enrichment'] = multipletests(cov_posctrl_ranks.P_Enrichment, \
                                                    alpha=0.25, method='fdr_bh')[1]
fdr = cov_posctrl_ranks.reindex(pos_ctrl_ranks.index).FDR_Enrichment
pos_ctrl_ranks['FDR_Enrichment'] = fdr
pos_ctrl_ranks.to_csv(outdir+'comparisons_validation_all-signatures.csv')

# Select validated comparisons
cmpr_sele = cov_posctrl_ranks.loc[(cov_posctrl_ranks['P_Enrichment'] < 0.05) & \
                                  (cov_posctrl_ranks['P_SpearmanR'] < 0.05) & \
                                  (cov_posctrl_ranks['FDR_Enrichment'] < 0.25) & \
                                  (cov_posctrl_ranks['SpearmanR_sRGES_EC50'] > 0.4)].index.to_list()
pos_ctrl_ranks.loc[cmpr_sele].to_csv(outdir+'comparisons_selected.csv')

# Plot the enrichment density and correlation score vs EC50
hnum = len(cmpr_sele) + len(cmpr_sele) % 2
FIG = plt.figure(figsize=(20, hnum * 2), dpi=300)
gs = gridspec.GridSpec(hnum, 4, height_ratios=[2 if i%2==0 else 1 for i in range(hnum)])
for i, comp in enumerate(cmpr_sele):
    # Edit title: remove internal code and add 'h' after time points
    title = comp[comp.index('GSE'):]
    tps = re.findall('_[0-9][0-9]*', title)
    for tp in tps: title = title.replace(tp, tp + 'h')
    # Find which sub panel
    if i%2 == 0: row, col = i, 0
    else: row, col = i - 1, 2
    # Draw density of positive drugs
    ax_dist = FIG.add_subplot(gs[row, col])
    p = pos_ctrl_ranks.loc[comp, 'P_Enrichment']
    sns.distplot(pos_ctrl_ranks.loc[comp, pos_ctrls], hist=False, \
                 kde_kws={'kernel': 'triw'}, ax=ax_dist, label='P = %.2E'%p)
    ax_dist.legend(fontsize=16)
    ax_dist.set_xlim(0, len(drug_lst))
    ax_dist.set_axis_off()
    ax_dist.set_title(title, fontsize=16)
    # Draw rank positions of positive drugs
    ax_bcd = FIG.add_subplot(gs[row + 1, col])
    ax_bcd.eventplot(pos_ctrl_ranks.loc[comp, pos_ctrls])
    ax_bcd.set_xlim(0, len(drug_lst))
    ax_bcd.set_ylim(0.5, 1)
    ax_bcd.set_xlabel('Positive Drugs Rank', fontsize=14)
    ax_bcd.get_yaxis().set_visible(False)
    for spine in ["left", "top", "right"]: ax_bcd.spines[spine].set_visible(False)
    # Draw correlation of srges and EC50
    ax_cor = FIG.add_subplot(gs[row:row + 2, col + 1])
    x = pos_ctrl_scores.loc[comp]
    y = pos_ctrl_ec50.loc[pos_ctrl_scores.columns, 'log10_EC50']
    corr = pos_ctrl_ranks.loc[comp, 'SpearmanR_sRGES_EC50']
    p = pos_ctrl_ranks.loc[comp, 'P_SpearmanR']
    sns.regplot(x, y, label='Spearman R = %.2f\nP = %.2E'%(corr, p), ax=ax_cor)
    ax_cor.set_xlabel('sRGES', fontsize=14)
    ax_cor.set_ylabel('log10 EC50', fontsize=14)
    ax_cor.set_title(title, fontsize=16)
    ax_cor.legend(fontsize=16)
FIG.tight_layout()
FIG.savefig(outdir + 'FigS4_Enrichment_EC50_sRGES.pdf', transparent=True)
plt.close()

# Plot the best signature result
comp = cmpr_sele[0]
title = comp[comp.index('GSE'):]
tps = re.findall('_[0-9][0-9]*', title)
for tp in tps: title = title.replace(tp, tp + 'h')
FIG = plt.figure(figsize=(10, 4), dpi=300)
gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1])
ax_dist = FIG.add_subplot(gs[0, 0])
p = pos_ctrl_ranks.loc[comp, 'P_Enrichment']
sns.distplot(pos_ctrl_ranks.loc[comp, pos_ctrls], hist=False, \
             kde_kws={'kernel': 'triw'}, ax=ax_dist, label='P = %.2E'%p)
ax_dist.legend(fontsize=16)
ax_dist.set_xlim(0, len(drug_lst))
ax_dist.set_axis_off()
ax_dist.set_title(title, fontsize=16)
ax_bcd = FIG.add_subplot(gs[1, 0])
ax_bcd.eventplot(pos_ctrl_ranks.loc[comp, pos_ctrls])
ax_bcd.set_xlim(0, len(drug_lst))
ax_bcd.set_ylim(0.5, 1)
ax_bcd.set_xlabel('Positive Drugs Rank', fontsize=14)
ax_bcd.get_yaxis().set_visible(False)
for spine in ["left", "top", "right"]: ax_bcd.spines[spine].set_visible(False)
ax_cor = FIG.add_subplot(gs[0:2, 1])
x = pos_ctrl_scores.loc[comp]
y = pos_ctrl_ec50.loc[pos_ctrl_scores.columns, 'log10_EC50']
corr = pos_ctrl_ranks.loc[comp, 'SpearmanR_sRGES_EC50']
p = pos_ctrl_ranks.loc[comp, 'P_SpearmanR']
sns.regplot(x, y, ax=ax_cor)
ax_cor.text(-0.2, -6.4, 'Spearman R = %.2f\nP = %.2E'%(corr, p), fontsize=16)
ax_cor.set_xlabel('sRGES', fontsize=14)
ax_cor.set_ylabel('log10 EC50', fontsize=14)
ax_cor.set_title(title, fontsize=16)
FIG.tight_layout()
FIG.savefig(outdir + 'Fig2A_Best_Enrichment_EC50_sRGES.pdf', transparent=True)
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
    if score_lst.sRGES.iloc[0] > -0.25: continue
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

# Combine activity table and toxicity/selectivity table
col_a = ['name', 'Median_Rank']
col_b = ['Median_Rank_MOCK'] + drug_info.columns.to_list()
com = pd.concat([conss_rank.loc[:, col_a], conss_rank_n.loc[:, col_b]], axis=1)
com.to_csv(outdir + 'consensus_rank_sRGES_drugs_CaseVSCtrl.csv', index=False)










