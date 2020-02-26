#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:30:38 2020

@author: jingxing
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'


inpdir = '../data/all_comparisons/'
outdir = '../data/consensus_rank/'
labelfile = outdir + 'comparisons_selected.csv'
lincs_sig = pd.read_hdf('/Users/jingxing/Documents/Data_Sources/LINCS_profiles/lincs_signatures_cmpd_landmark_highQ.hdf5')
meta = pd.read_csv('/Users/jingxing/Documents/Data_Sources/LINCS_profiles/signature_meta.csv', index_col='id')

pos_drug_ranks = pd.read_csv(labelfile, index_col=0)
tmp = pos_drug_ranks.iloc[:, 1:-3].T
tmp['med'] = tmp.median(axis=1)
pos_ctrls = tmp.sort_values(by='med').index.to_list()
models = pos_drug_ranks.index.to_list()
sele_drugs = ['bortezomib', 'puromycin', 'methotrexate', 'methylene-blue', 'tyloxapol', \
              'nisoldipine', 'nvp-bez235', 'fluvastatin', 'alvocidib', 'dasatinib', \
              'orteronel', 'bezafibrate', 'atovaquone']
colorp = 'coolwarm'
use_median = False

# Summarize selected dz signatures
gene_all = set()
for m in models:
    filename = inpdir + m + '/dz_sig_used.csv'
    dzus = pd.read_csv(filename)
    gene_all = gene_all.union(dzus.Symbol)
gene_all = sorted(gene_all)
sigmtx = []
for m in models:
    filename = inpdir + m + '/dz_sig_used.csv'
    dzsig = pd.read_csv(filename, index_col='Symbol')
    sigmtx.append([dzsig.loc[g, 'log2FoldChange'] if g in dzsig.index else 0.0 for g in gene_all])
sigmtx = np.array(sigmtx).T
sigmtx = pd.DataFrame(data=sigmtx, index=gene_all, columns=models)
# Calculate consensus signature
sigmtx['Num_Pos'] = sigmtx.where(sigmtx > 0).count(axis=1)
sigmtx['Num_Neg'] = sigmtx.where(sigmtx < 0).count(axis=1)
sigmtx['Num_NonZero'] = sigmtx.where(sigmtx != 0).count(axis=1)
sigmtx['Q25'] = sigmtx.quantile(0.25, axis=1)
sigmtx['Q75'] = sigmtx.quantile(0.75, axis=1)
sigmtx_cons = sigmtx.loc[(sigmtx.Q75 > 2) | (sigmtx.Q25 < -1)]
sigmtx_cons['Infection'] = [a if a < 0 else b * 0.75 for (a,b) in \
                            zip(sigmtx_cons.Q25, sigmtx_cons.Q75)]

def funs(score_lst):
    x = np.quantile(score_lst, 0.25)
    y = np.quantile(score_lst, 0.75)
    if x > -1 and y < 1: return np.median(score_lst)
    if x < -y: return np.median([i for i in score_lst if i <= x])
    else: return np.median([i for i in score_lst if i >= y])
    

# Collect drug profiles
meta.index = meta.index.astype('str')
meta = meta.loc[meta.index.isin(lincs_sig.columns)]
meta.pert_iname = meta.pert_iname.str.lower()
gene_sele = sorted(sigmtx_cons.index)
lincs_sig = lincs_sig.loc[gene_sele]


# Plot selected drugs
sub_meta = meta.loc[(meta.pert_iname.isin(sele_drugs)) & (meta.pert_dose == 10)]
sele_drugs = [d for d in sele_drugs if d in sub_meta.pert_iname.values]
sele_drug_profiles = []
for drug in sele_drugs:
    prf = lincs_sig.loc[:, lincs_sig.columns.isin(sub_meta.loc[sub_meta.pert_iname == drug].index)]
    if use_median:
        m_prf = prf.median(axis=1)
        sele_drug_profiles.append(m_prf.to_list())
    else:
        m_prf = [funs(t) for k,t in prf.iterrows()]
        sele_drug_profiles.append(m_prf)
    
sele_drug_profiles = np.array(sele_drug_profiles).T
sele_drug_profiles = pd.DataFrame(data=sele_drug_profiles, columns=sele_drugs, index=lincs_sig.index)
dz_sele_drug = pd.concat([sigmtx_cons['Infection'], sele_drug_profiles], axis=1)
dz_sele_drug = dz_sele_drug.sort_values(by='Infection')
dz_sele_drug.index.name = 'Dysregulated Genes'
FIG = plt.figure(figsize=(14, 20), dpi=300)
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 12])
ax0 = FIG.add_subplot(gs[0, 0])
sns.heatmap(dz_sele_drug.iloc[:, [0]], cmap=colorp, cbar=False, ax=ax0)
ax0.set_xticklabels(['Infection'], rotation=45, fontsize=24)
ax0.set_yticklabels(dz_sele_drug.index, fontsize=6)
ax0.set_ylabel(dz_sele_drug.index.name, fontsize=24)
ax1 = FIG.add_subplot(gs[0, 1])
sns.heatmap(dz_sele_drug.iloc[:, 1:], cmap=colorp, cbar=False, ax=ax1)
ax1.get_yaxis().set_visible(False)
ax1.set_xticklabels(sele_drug_profiles.columns, rotation=45, fontsize=24)
FIG.tight_layout()
FIG.savefig(outdir + 'Fig2C_SigRev_SelectedDrugs.pdf', transparent=True)
plt.close()

FIG = sns.clustermap(dz_sele_drug, cmap=colorp, figsize=(10, 30), dendrogram_ratio=(0.1, 0.05))
FIG.savefig(outdir + 'SigRevCluster_SelectedDrugs.pdf')
plt.close()







