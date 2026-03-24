#!/usr/bin/env python3


import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# TOP MUTATED GENES
print("HNSC TOP MUTATED GENES")
print("=" * 50)


# Load cleaned data
mutations_clean = pd.read_csv('../data/clean_mutations.tsv', sep='\t')
clinical_clean = pd.read_csv('../data/clean_clinical.tsv', sep='\t')


print(f"Loaded: {mutations_clean.shape} mutations, {clinical_clean.shape} clinical")
print("Sample of data:")
print(mutations_clean.head(3))


# Exlude hypermutators
non_hyper = mutations_clean[mutations_clean['hypermutator'] != 1]
print(f"\nSamples before hypermutator filter: {mutations_clean['SAMPLE_ID'].nunique()}")
print(f"Samples after: {non_hyper['SAMPLE_ID'].nunique()}")

total_samples = non_hyper['SAMPLE_ID'].nunique()
print(f"Denominator: {total_samples} non-hypermutated primary samples")



# Gene frequencies (unique samples per gene)
gene_freq = (non_hyper.groupby('Hugo_Symbol')['SAMPLE_ID']
                      .nunique()
                      .transform(lambda x: x / total_samples * 100)
                      .sort_values(ascending=False))   

print("\nTOP 30 GENE FREQUENCIES:")
top30 = gene_freq.head(30)
print(top30.round(1))



# Validation check
drivers = ['TP53', 'CDKN2A', 'FAT1', 'NOTCH1', 'PIK3CA']
driver_ranks = [top30.index.get_loc(gene) + 1 for gene in drivers if gene in top30]
print(f"\nVALIDATION - Drivers in top 15? {all(r <= 15 for r in driver_ranks)}")
print("Driver ranks:", dict(zip(drivers, driver_ranks)))


# FIGURE: Top 25 mutated genes
top25_genes = top30.head(25).index
top25_data = top30.head(25)



# Drivers vs passengers (manual curation from Day 3)
known_drivers = ['TP53', 'CDKN2A', 'FAT1', 'NOTCH1', 'PIK3CA', 'KMT2D', 'CASP8', 'NSD1']
colors = ['darkred' if gene in known_drivers else 'lightgreen' for gene in top25_genes]

plt.figure(figsize=(10, 12))
bars = top25_data.plot(kind='barh', color=colors, edgecolor='black', linewidth=0.5)
plt.title(f'Top 25 Mutated Genes in HNSC\nNon-silent mutations, hypermutators excluded\n(n={total_samples} samples)', fontsize=14, fontweight='bold', pad=20)
plt.xlabel('% Non-hypermutated Primary Samples Mutated', fontsize=12)
plt.xlim(0, 80)
plt.axvline(x=15, color='red', linestyle='--', alpha=0.7, label='15% threshold')
legend_elements = [
    Patch(facecolor='darkred', edgecolor='black', label='Known drivers'),
    Patch(facecolor='lightgreen', edgecolor='black', label='Passenger/Large genes'),
    Line2D([0], [0], color='red', linestyle='--', label='15% threshold')
]
plt.legend(handles=legend_elements, loc='center right', fontsize=10)
plt.tight_layout()


# Add % labels
for i, (gene, pct) in enumerate(zip(top25_genes, top25_data)):
    plt.text(pct + 0.5, i, f'{pct:.1f}%', va='center', fontweight='bold')

plt.savefig('../figures/mutation_frequency.png', dpi=300, bbox_inches='tight')
plt.show()


print("Figure saved: ../figures/mutation_frequency.png")