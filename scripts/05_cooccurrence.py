#!/usr/bin/env python3


"""
HNSC Co-occurrence Heatmap (Production)
Pairwise Fisher's exact + BH correction
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from itertools import combinations
from statsmodels.stats.multitest import multipletests
import os


# Paths
os.makedirs('../figures', exist_ok=True)


print("HNSC CO-OCCURENCE HEATMAP")
print("@%"* 50)



# Load mutations
mutations = pd.read_csv('../data/clean_mutations.tsv', sep='\t')
non_hyper_mut = mutations[mutations['hypermutator'] !=1]
all_samples = non_hyper_mut['SAMPLE_ID'].unique()
total_samples = len(all_samples)
print(f"Samples: {total_samples}")


# Top 15 genes
top_genes = non_hyper_mut.groupby('Hugo_Symbol')['SAMPLE_ID'].nunique().sort_values(ascending=False).head(15).index
print(f"Genes: {len(top_genes)}")


# Binary matrix (504 x 15 genes)
binary_data = non_hyper_mut[non_hyper_mut['Hugo_Symbol'].isin(top_genes)]
binary_pivot = binary_data.pivot_table(
    index='SAMPLE_ID', columns='Hugo_Symbol', values='Variant_Classification', aggfunc='size', fill_value=0
).gt(0).astype(int)
binary_pivot = binary_pivot.reindex(all_samples, fill_value=0)
print(f"Binary matrix: {binary_pivot.shape}")


# All pairwise Fisher's exact tests
print("\nComputing pairwise tests...")
results = []
for g1, g2 in combinations(top_genes, 2):
    col1 = binary_pivot[g1]
    col2 = binary_pivot[g2]


    a = ((col1 == 1) & (col2 == 1)).sum()
    b = ((col1 == 1) & (col2 == 0)).sum()
    c = ((col1 == 0) & (col2 == 1)).sum()
    d = ((col1 == 0) & (col2 == 0)).sum()
    
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
    
    results.append({
        'gene1': g1, 'gene2': g2, 
        'odds_ratio': odds_ratio, 'p_value': p_value
    })


results_df = pd.DataFrame(results)


# BH correction
rejected, q_values, _, _ = multipletests(results_df['p_value'], method='fdr_bh', alpha=0.05)
results_df['q_value'] = q_values
results_df['significant'] = rejected


print(f"Significant pairs (q<0.05): {results_df['significant'].sum()}/{len(results_df)}")



# HEATMAP MATRIX
print("\nBuilding heatmap")
genes = top_genes.tolist()
n = len(genes)
plot_matrix = pd.DataFrame(0.0, index=genes, columns=genes)



for _, row in results_df.iterrows():
    score = -np.log10(row['q_value'] + 1e-10)       # -log10(q)
    if row['odds_ratio'] < 1:
        score = -score          # Negative = mutual exclusivity
    plot_matrix.loc[row['gene1'], row['gene2']] = score
    plot_matrix.loc[row['gene2'], row['gene1']] = score


# Mask upper triangle
mask = np.triu(np.ones_like(plot_matrix, dtype=bool))
plot_matrix_display = plot_matrix.where(~mask, np.nan)

plt.figure(figsize=(12, 10))
sns.heatmap(plot_matrix_display, cmap='RdBu_r', center=0, annot=False, square=True, cbar_kws={'label': '-log10(q) x direction'}, xticklabels=True, yticklabels=True, mask=mask)


plt.title(f'HNSC Gene Co-occurrence/Mutual Exclusivity\n' f'Top 15 genes, {total_samples} samples, BH-corrected (q<0.05)',
fontsize = 14, fontweight='bold', pad=20)
plt.xlabel('Gene 2', fontsize=12)
plt.ylabel('Gene 1', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)


# Annotate signification pairs
# Annotate significant pairs (lower triangle: row > col)
for i, (g1, g2) in enumerate(combinations(genes, 2)):
    if results_df.iloc[i]['significant']:
        i_idx = genes.index(g1)
        j_idx = genes.index(g2)
        r, c = max(i_idx, j_idx), min(i_idx, j_idx)
        plt.text(c + 0.5, r + 0.5, '*', fontsize=14, ha='center', va='center',
                 color='black', fontweight='bold')

plt.tight_layout()
plt.savefig('../figures/cooccurrence.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()


print("\n SAVED: figures/cooccurrence.png")
print("\nSUMMARY:")
print(results_df[results_df['significant']][['gene1', 'gene2', 'odds_ratio', 'q_value']].round(3).to_string())