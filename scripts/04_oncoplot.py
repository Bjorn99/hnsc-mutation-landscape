#!/usr/bin/env python3

"""
From notebooks/oncoplot_dev.ipynb matrix
"""


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import numpy as np 
import os


# Paths
os.makedirs('../figures', exist_ok=True)

print("HNSC ONCOPLOT - PRODUCTION")
print("%" * 50)


# Load Matrix
matrix = pd.read_csv('../data/oncoplot_matrix.csv', index_col=0)
print(f"Matrix: {matrix.shape}")



# Mutation type colors
colors = {
    'Missense_Mutation': '#26A269',      # Green
    'Nonsense_Mutation': '#2E2E2E',      # Black
    'Frame_Shift_Del': '#2A7FFF',        # Blue
    'Frame_Shift_Ins': '#E66100',        # Orange
    'Splice_Site': '#E6C700',            # Yellow
    'In_Frame_Del': '#D4726A',           # Pink
    'In_Frame_Ins': '#8B4513',           # Brown
    'Translation_Start_Site': '#CC79A7', # Pink
    'Nonstop_Mutation': '#F0E442'        # Yellow
}


# Create numeric matrix for imshow
unique_types = sorted(matrix.stack().dropna().unique())
type_to_num = {t: i+1 for i, t in enumerate(unique_types)}  #1=mutated, 0=WT


# Custom colormap (WT=white, mutations=correct colors)
ordered_types = [''] + list(type_to_num.keys())     # '' for WT=0
cmap_colors = ['white'] + [colors.get(t, 'lightgray') for t in type_to_num.keys()]
custom_cmap = ListedColormap(cmap_colors)


num_matrix = matrix.map(lambda x: type_to_num.get(x, 0)).values

print(f"Unique mutation types: {len(unique_types)}")
print("Colormap ready:, WT=white")


# PLOT
fig, (ax, ax_freq) = plt.subplots(2, 1, figsize=(22, 12), gridspec_kw={'height_ratios': [20, 1]}, sharex=True)

# Oncoplot heatmap
im = ax.imshow(num_matrix, aspect='auto', cmap=custom_cmap, interpolation='none', vmin=0, vmax=len(unique_types)+1)
ax.set_yticks(np.arange(len(matrix.index)))
ax.set_yticklabels(matrix.index, fontsize=10, fontweight='bold')


# Sample mutation counts (top bar)
sample_mut_counts = matrix.notna().sum(axis=0)
ax_freq.bar(range(len(matrix.columns)), sample_mut_counts.values, color='gray', alpha=0.6, width=1)
ax_freq.set_ylabel('Mutations/Sample', fontsize=11)
ax_freq.set_ylim(0, sample_mut_counts.max() * 1.1)


# Formatting
ax.set_title('HNSC Oncoplot: Top 20 Mutated Genes\n'
f'Non-silent mutations, hypermutators excluded (n = {matrix.shape[1]} of 504 samples)', fontsize=16, fontweight='bold', pad=20)

ax.set_xlabel('Samples (sorted by TP53 status)', fontsize=12)


# Legend matches colors exactly
patches = [mpatches.Patch(color=colors[t], label=t) for t in unique_types[:8] if t in colors]
patches.append(mpatches.Patch(color='white', label='Wildtype'))
ax.legend(handles=patches, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)

plt.tight_layout()
plt.savefig('../figures/oncoplot_fixed.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()



print(f"\n SAVED: figures/oncoplot_fixed.png")
print(f" TP53 cluster visible (left: {matrix.loc['TP53'].notna().mean()*100:.1f}%)")