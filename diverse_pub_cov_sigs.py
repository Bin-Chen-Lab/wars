#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 20:04:37 2022

@author: jingxing
"""

import pandas as pd
import os
from scipy.stats import ranksums, spearmanr, fisher_exact
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
import seaborn as sns


inpdir = '../data/published_covid_signatures/'
outdir = '../data/published_covid_signatures/SUMMARY/'
pub_cov_cpd = pd.read_csv(outdir + 'pub_cov_cmpds.csv', index_col=0)
drug_lib = pd.read_csv('/Users/jingxing/Documents/OneDrive - Michigan State University/Data_Sources/LINCS_cmpds/repurposing_drugs_20180907_smiles.csv', index_col='pert_iname')

rges_all, rges_rank_all, DE_all = [], [], []
for fld in os.listdir(inpdir):
    if fld.startswith('.') or fld == 'SUMMARY': continue
    print(fld)
    tmp = pd.read_csv(inpdir + fld + '/sRGES.csv', index_col='pert_iname')
    tmp = tmp.loc[(tmp.n > 1) & (tmp.index.isin(drug_lib.index)), 'sRGES']
    tmp.name = fld
    rges_all.append(tmp)
    tmp = tmp.rank()
    rges_rank_all.append(tmp)
    tmp = pd.read_csv(inpdir + fld + '/DZSIG__%s.csv'%fld, index_col='Symbol')
    tmp = tmp.groupby(by=tmp.index).mean()
    tmp[fld] = [1 if lfc > 0 else -1 for lfc in tmp['log2FoldChange']]
    DE_all.append(tmp[fld])
rges_all = pd.concat(rges_all, axis=1)
rges_all = rges_all.loc[:, sorted(rges_all.columns)]
rges_rank_all = pd.concat(rges_rank_all, axis=1)
rges_rank_all = rges_rank_all.loc[:, sorted(rges_rank_all.columns)]
DE_all = pd.concat(DE_all, axis=1)
DE_all = DE_all.loc[:, sorted(DE_all.columns)]

rges_ec50_all = pub_cov_cpd[pub_cov_cpd.index.isin(rges_all.index)]
rges_ec50_all = pd.concat([rges_ec50_all, rges_all.reindex(rges_ec50_all.index)], axis=1)
actives = set(rges_ec50_all[rges_ec50_all['log10_EC50_v2'] < -5].index)

m, n, libsize = 3, 6, rges_rank_all.shape[0]
FIG, axes = plt.subplots(m, n, figsize=(n * 2.16, m * 1.89))
for num, dzsig in enumerate(rges_all.columns):
    ay, ax = int(num / n), num % n
    cor, p = spearmanr(rges_ec50_all['log10_EC50_v2'], rges_ec50_all[dzsig])
    CLR = 'black'
    if p < 0.01 and cor > 0: CLR = 'royalblue'
    if p < 0.01 and cor < 0: CLR = 'grey'
    sns.regplot(x=dzsig, y='log10_EC50_v2', data=rges_ec50_all, ax=axes[ay, ax], color=CLR)
    if ay == m - 1:
        axes[ay, ax].set_xlabel('Compound sRGES')
    else:
        axes[ay, ax].set_xlabel(None)
    if ax == 0:
        axes[ay, ax].set_ylabel('Antiviral log10 EC50')
    else:
        axes[ay, ax].set_ylabel(None)
    axes[ay, ax].set_title('%s\nCor = %.2f'%(dzsig, cor))
plt.tight_layout()
FIG.savefig(outdir + 'Cor_sRGES_EC50_pub_cov_sigs.pdf', transparent=True)

FIG, axes = plt.subplots(m * 2, n, figsize=(n * 3, m * 3), \
                         gridspec_kw={'height_ratios': [3 if i%2==0 else 1 for i in range(m * 2)]})
for num, dzsig in enumerate(rges_all.columns):
    ay0, ax = int(num / n) * 2, num % n
    ay1 = ay0 + 1
    pos_ranks = rges_rank_all.loc[rges_rank_all.index.isin(actives), dzsig]
    _, p = ranksums(pos_ranks, rges_rank_all.loc[~rges_rank_all.index.isin(rges_ec50_all.index), dzsig])
    p = p * 0.5 if _ < 0 else 1 - p * 0.5
    CLR = 'black'
    if p < 0.01: CLR = 'royalblue'
    if p > 0.99: CLR = 'grey'
    ptxt = 'P = %.2f'%p if p > 0.01 else 'P = %.2E'%p
    sns.distplot(pos_ranks, hist=False, kde_kws={'kernel': 'triw'}, ax=axes[ay0, ax], color=CLR)
    axes[ay0, ax].set_xlim(0, libsize)
    axes[ay0, ax].set_axis_off()
    axes[ay0, ax].set_xlabel(None)
    axes[ay0, ax].set_title(dzsig + '\n' + ptxt)
    axes[ay1, ax].eventplot(pos_ranks, color=CLR)
    axes[ay1, ax].set_xlim(0, libsize)
    axes[ay1, ax].set_ylim(0.5, 1)
    if ay1 == m * 2 - 1:
        axes[ay1, ax].set_xlabel('Rank of SARS-CoV-2 Inhibitors')
    else:
        axes[ay1, ax].set_xlabel(None)
    axes[ay1, ax].get_yaxis().set_visible(False)
    for spine in ["left", "top", "right"]: axes[ay1, ax].spines[spine].set_visible(False)
plt.tight_layout()
FIG.savefig(outdir + 'Enrichment_sRGES_actives_pub_cov_sigs.pdf', transparent=True)

t = DE_all.count(axis=1)
viz = DE_all.loc[t >= 3]
viz = viz.fillna(0).T
FIG = sns.clustermap(viz, center=0, cmap='coolwarm', figsize=(10, 6), xticklabels=False, dendrogram_ratio=(0.1, 0.1))
FIG.savefig(outdir + 'Cluster_pub_cov_sigs.pdf', transparent=True)

gene_set_file = '/Users/jingxing/Library/CloudStorage/OneDrive-MichiganStateUniversity/Data_Sources/Gene_sets/MSigDB_Hallmark_2020.txt'
gene_sets = []
with open(gene_set_file) as f:
    for line in f:
        key, l = line.split('\t\t')
        gene_sets.append([key, set(l.split('\t'))])
gene_sets = sorted(gene_sets)
univ = set.union(*[v for k, v in gene_sets])
N = len(univ)
pw = []
for dzsig in DE_all.columns:
    dz_up = set(DE_all.loc[DE_all[dzsig] > 0].index)
    dz_down = set(DE_all.loc[DE_all[dzsig] < 0].index)
    n_up, n_down = len(dz_up), len(dz_down)
    row = [dzsig]
    for key, gs in gene_sets:
        k_up = len(set.intersection(dz_up, gs))
        K = len(gs)
        _, p_up = fisher_exact([[k_up, K], [n_up, N]], alternative='greater')
        k_down = len(set.intersection(dz_down, gs))
        _, p_down = fisher_exact([[k_down, K], [n_down, N]], alternative='greater')
        logp = -np.log10(p_up) if p_up < p_down else np.log10(p_down)
        row.append(logp)
    pw.append(row)
pw = pd.DataFrame(data=pw, columns=['COVID Signature'] + [k for k, v in gene_sets])
pw.index = pw.iloc[:, 0]
pw = pw.iloc[:, 1:]
count_col = pw.where(pw.abs() >= -np.log10(0.05)).count()
viz = pw.loc[:, count_col >= 1].T

FIG = sns.clustermap(viz, center=0, vmin=-2, vmax=2, cmap='coolwarm', figsize=(6.4, 4.8), dendrogram_ratio=(0.2, 0.1), cbar_pos=(0.02, 0.4, 0.05, 0.18))
ax = FIG.ax_heatmap
ax.set_xticklabels(ax.get_xticklabels(),rotation = 45, ha='right')
#ax.set_ylabel('MSigDB Hallmark')
ax.set_xlabel('Published COVID-19 Signatures', fontsize=12)
FIG.savefig(outdir + 'MSigDB_enrichment_pub_dzsigs_cluster.pdf', transparent=True)


