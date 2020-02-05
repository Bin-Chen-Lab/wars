#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 09:55:26 2019

@author: jingxing
"""

import json
from itertools import combinations
import matplotlib.pyplot as plt
import os
import pandas as pd

# Summarize a few disease signatures and generate a consensus disease signature
inpdir = '/Users/jingxing/Documents/CoV/WARS/data/comparisons/'
dz_flds = pd.read_csv('/Users/jingxing/Documents/CoV/WARS/data/consensus_rank/comparisons_selected.csv').model.to_list()
#dz_flds = [fld for fld in os.listdir(inpdir) if not fld.startswith('.')]
dz_sig_files = [inpdir + fld + '/dz_sig_used.csv' for fld in dz_flds]
gene_all = set()
for fname in dz_sig_files:
    gene_all = gene_all.union(pd.read_csv(fname).Symbol)
gene_all = sorted(gene_all)
dz_sig_dict = dict([(g, []) for g in gene_all])
for fname in dz_sig_files:
    dz_sig = pd.read_csv(fname, index_col='Symbol')
    for gene in dz_sig_dict:
        foldchg = dz_sig.loc[gene, 'value'] if gene in dz_sig.index else 0.0
        dz_sig_dict[gene].append(foldchg)
dz_sig_c = [dz_sig_dict[g] for g in gene_all]
dz_sig_c = pd.DataFrame(data=dz_sig_c, index=gene_all, columns=dz_flds)
dz_sig_c['Ave_CoV_Signature'] = dz_sig_c.mean(axis=1)
    

# Load PharmaocoDepMap
pdm = json.load(open('../data/raw/PharmacoDepMap_all_d1t0_200117.json'))
# Lung cell lines
cell_lines = ['A549', 'HCC515', 'DV90', 'H1299', 'HCC15', 'NCIH596', 'NCIH2073', 'SKLU1']
outdir = '../data/combination/'

result = dict([(dz, {1: [], 2:[]}) for dz in dz_sig_c.columns])

for dz in dz_sig_c.columns:
    # Summarize disease signature reversal dependency
    dz_sig = dict([(g, v) for g,v in zip(dz_sig_c.index, dz_sig_c[dz]) if v != 0])
    dzgene2tar = dict()
    for gene, fldchg in dz_sig.items():
        dr = 'Up' if fldchg < 0 else 'Down'
        tlst = []
        for cell_line in cell_lines:
            if cell_line in pdm[dr] and gene in pdm[dr][cell_line]:
                tlst += pdm[dr][cell_line][gene]
        if gene in pdm[dr]['Universal']:
            tlst += pdm[dr]['Universal'][gene]
        tlst = list(set(tlst))
        dzgene2tar[gene] = tlst
    tar2dzgene = dict()
    for gene, tlst in dzgene2tar.items():
        for t in tlst:
            if t in tar2dzgene:
                tar2dzgene[t].append(gene)
            else:
                tar2dzgene[t] = [gene]
    
    # Search for combinations
    single_target_effects = []
    for t in tar2dzgene:
        dz_genes = tar2dzgene[t]
        s = sum([abs(dz_sig[rq]) for rq in dz_genes])
        single_target_effects.append([t, s, len(dz_genes), '|'.join(dz_genes)])
        if len(single_target_effects) > 10:
            single_target_effects = sorted(single_target_effects, key=lambda x:x[1], reverse=True)[:10]
    result[dz][1] = single_target_effects
    double_targets_effects = []
    for targets in combinations(tar2dzgene, 2):
        dz_genes = set.union(*[set(tar2dzgene[t]) for t in targets])
        s = sum([abs(dz_sig[rq]) for rq in dz_genes])
        double_targets_effects.append(['--'.join(targets), s, len(dz_genes), '|'.join(dz_genes)])
        if len(double_targets_effects) > 10:
            double_targets_effects = sorted(double_targets_effects, key=lambda x:x[1], reverse=True)[:10]
    result[dz][2] = double_targets_effects

    # Output top 10 single and double targets into csv
    table = single_target_effects + double_targets_effects
    table = sorted(table, key=lambda x:x[1], reverse=True)
    with open(outdir + 'ReversalTargets_%s.csv'%dz, 'w') as f:
        f.write('Targets,Score,Num_Rev_Genes,Rev_Genes\n')
        for t, s, num, gl in table:
            f.write('%s,%.2f,%s,%s\n'%(t, s, num, gl))
    
    best = max(single_target_effects + double_targets_effects, key=lambda x:x[1])
    
    # Waterfall plot of the highest reversal score combination
    FIG = plt.figure(figsize=(12,7), dpi=300)
    dz_gene_wt = sorted([(k, v) for k,v in dz_sig.items()], key=lambda x:x[1])
    X0 = [i[0] for i in dz_gene_wt]
    Y0 = [i[1] for i in dz_gene_wt]
    plt.bar(x = X0, height=Y0, label=dz, color='grey')
    rev_gene_wt = sorted([(k, dz_sig[k]) for k in best[-1].split('|')])
    X1 = [i[0] for i in rev_gene_wt]
    Y1 = [i[1] for i in rev_gene_wt]
    plt.bar(x = X1, height=Y1, label=best[0]+' reversed genes', color='red')
    plt.legend(fontsize=14)
    plt.xticks(rotation=90, fontsize=10)
    plt.ylabel('Log2 Fold Change (Disease vs Control)', fontsize=12)
    FIG.tight_layout()
    FIG.savefig(outdir + '%s_%s.pdf'%(dz, best[0]))
    plt.close()











    