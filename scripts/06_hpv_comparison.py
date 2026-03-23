#!/usr/bin/env python3

"""
HPV+ vs HPV- Mutation Comparison
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import os

os.makedirs('../figures', exist_ok=True)

print("HPV+ vs HPV- MUTATION LANDSCAPE")
print("%" * 50)


# Load + filter data
mutations = pd.read_csv('../data/clean_mutations.tsv', sep='\t')
clinical = pd.read_csv('../data/clean_clinical.tsv', sep='\t')


non_hyper_mut = mutations[mutations['hypermutator'] != 1]
non_hyper_clin = clinical[clinical['hypermutator'] != 1]


# HPV groups (samples WITH mutations)
hpv_pos = non_hyper_clin[non_hyper_clin['HPV_clean'] == 'Positive']['SAMPLE_ID']
hpv_neg = non_hyper_clin[non_hyper_clin['HPV_clean'] == 'Negative']['SAMPLE_ID']

hpv_pos_mut = hpv_pos[hpv_pos.isin(non_hyper_mut['SAMPLE_ID'])]
hpv_neg_mut = hpv_neg[hpv_neg.isin(non_hyper_mut['SAMPLE_ID'])]


n_pos, n_neg = len(hpv_pos_mut), len(hpv_neg_mut)
print(f"HPV+ (n={n_pos}), HPV- (n={n_neg})")


# TOP 15 GENES
top_genes = non_hyper_mut.groupby('Hugo_Symbol')['SAMPLE_ID'].nunique().sort_values(ascending=False).head(15).index.tolist()


# Frequencies + Fisher's tests
results = []
for gene in top_genes:
    pos_mut = non_hyper_mut[(non_hyper_mut['SAMPLE_ID'].isin(hpv_pos_mut)) & (non_hyper_mut['Hugo_Symbol'] == gene)]['SAMPLE_ID'].nunique()
    neg_mut = non_hyper_mut[(non_hyper_mut['SAMPLE_ID'].isin(hpv_neg_mut)) & (non_hyper_mut['Hugo_Symbol'] == gene)]['SAMPLE_ID'].nunique()
    freq_pos, freq_neg = pos_mut/n_pos*100, neg_mut/n_neg*100

    # Fisher's test
    a, b, c, d = pos_mut, n_pos-pos_mut, neg_mut, n_neg-neg_mut
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]])

    results.append({
        'gene': gene, 'hpv_pos_freq': freq_pos, 'hpv_neg_freq': freq_neg, 'odds_ratio': odds_ratio, 'p_value': p_value
    })

results_df = pd.DataFrame(results)


# BH correction
rejected, q_values, _, _ = multipletests(results_df['p_value'], method='fdr_bh', alpha=0.05)
results_df['q_value'] = q_values
results_df['significant'] = rejected

print(f"Significant genes: {results_df['significant'].sum()}/15")
print(results_df[results_df['significant']][['gene', 'hpv_pos_freq', 'hpv_neg_freq']].round(1))


# PLOT: Horizontal grouped bars
results_df['freq_diff'] = results_df['hpv_pos_freq'] - results_df['hpv_neg_freq']
results_df = results_df.sort_values('freq_diff').reset_index(drop=True)

fig, ax = plt.subplots(figsize=(10, 12))

x_pos = np.arange(len(results_df))
width = 0.35


bars1 = ax.barh(x_pos - width/2, results_df['hpv_pos_freq'], width, label=f'HPV+ (n={n_pos})', color='#1f77b4', alpha=0.8)
bars2 = ax.barh(x_pos + width/2, results_df['hpv_neg_freq'], width, label=f'HPV- (n={n_neg})', color='#ff7f0e', alpha=0.8)


# Gene labels
ax.set_title(f'HPV+ vs HPV- Mutation Frequencies in HNSC\n'
    f'Non-silent mutations, BH-corrected (* q<0.05)',
    fontsize=14, fontweight='bold', pad=20)
ax.set_xlabel('% Samples Mutated', fontsize=12)
ax.set_yticks(x_pos)
ax.set_yticklabels(results_df['gene'], fontweight='bold', fontsize=11)


# Asterisks for significant
for i, row in results_df.iterrows():
    if row['significant']:
        ax.text(90, i, '*', fontsize=18, fontweight='bold', ha='center', va='center', color='black')

ax.legend(loc='upper right')
ax.grid(axis='x', alpha=0.3)
ax.set_xlim(0, 95)


plt.tight_layout()
plt.savefig('../figures/hpv_comparison.png', dpi=300, bbox_inches='tight')
plt.show()


print("\n SAVED: figures/hpv_comparison.png")