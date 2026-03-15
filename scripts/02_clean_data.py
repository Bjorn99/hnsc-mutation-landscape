#!/usr/bin/env python3


import pandas as pd
import os


DATA_DIR = "../data/hnsc_tcga"
files = {
        'mutations': os.path.join(DATA_DIR, "data_mutations.txt"),
        'patient': os.path.join(DATA_DIR, "data_clinical_patient.txt"),
        'sample': os.path.join(DATA_DIR, "data_clinical_sample.txt")
}

# Load all files with pandas defaults
mutations = pd.read_csv(files['mutations'], sep='\t',
                        low_memory=False)
patient_data = pd.read_csv(files['patient'], sep='\t', comment='#',
                           low_memory=False)
sample_data = pd.read_csv(files['sample'], sep='\t', comment='#',
                          low_memory=False)


print(f"Raw: mutations {mutations.shape[0]:,}, samples {sample_data.shape}, patients {patient_data.shape}")

# CLEAN MISSING VALUES
print("\n Cleaning [Not Available]...")
na_vals = ['[Not Available]', '[Not Evaluated]', '[Not Applicable]', '[Discrepancy]']
mutations = mutations.replace(na_vals, pd.NA)
sample_data = sample_data.replace(na_vals, pd.NA)
patient_data = patient_data.replace(na_vals, pd.NA)


# NUMERIC CONVERSION
print("\n Numeric conversion...")
patient_data['OS_MONTHS'] = pd.to_numeric(patient_data['OS_MONTHS'], errors='coerce')
patient_data['DFS_MONTHS'] = pd.to_numeric(patient_data['DFS_MONTHS'], errors='coerce')
patient_data['AGE'] = pd.to_numeric(patient_data['AGE'], errors='coerce')
print("OS_MONTHS dtype:", patient_data['OS_MONTHS'].dtype)


# DUPLICATE SAMPLES
print("\n Duplicate samples check...")
dup_patients = sample_data[sample_data.duplicated('PATIENT_ID', keep=False)]
print("Patients with multiple samples:")
print(dup_patients[['PATIENT_ID', 'SAMPLE_ID', 'SAMPLE_TYPE']].drop_duplicates())


# Keep PRIMARY TUMOR (-01 suffix only)
sample_data['suffix'] = sample_data['SAMPLE_ID'].str[-2:]
primary_samples = sample_data[sample_data['suffix'] == '01'].copy()
print(f"Samples: {len(sample_data)} → Primary: {len(primary_samples)}")


# FILTER MUTATION TYPES
print("\n Non-silent mutation only...")
coding_types = [
        'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
        'In_Frame_Del', 'In_Frame_Ins', 'Translation_Start_Site', 'Splice_Site', 'Nonstop_Mutation'
]
mutations_clean = mutations[mutations['Variant_Classification'].isin(coding_types)].copy()
print(f"Mutations: {len(mutations):,} -> Coding: {len(mutations_clean):,}({100*(1-len(mutations_clean)/len(mutations)):.1f}% filtered)")


# SELECT COLUMNS
print("\n Select key columns...")
keep_cols = [
        'Hugo_Symbol', 'Entrez_Gene_Id', 'Variant_Classification', 'Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'End_Position'
]
if 'HGVSp_Short' in mutations_clean.columns:
        keep_cols.append('HGVSp_Short')
mutations_final = mutations_clean[keep_cols].copy()
print(f"Columns: {len(keep_cols)} selected")

# KEEP ONLY PRIMARY TUMOR MUTATIONS
mutations_final = mutations_final[
    mutations_final['Tumor_Sample_Barcode'].isin(primary_samples['SAMPLE_ID'])
].copy()
print(f"After primary filter: {len(mutations_final):,} mutations")


# HYPERMUTATORS
print("\n Hypermutators (filtered data)...")
sample_counts = mutations_final['Tumor_Sample_Barcode'].value_counts()
threshold = sample_counts.mean() + 3 * sample_counts.std()
hypermutators = sample_counts[sample_counts > threshold].index

primary_samples['hypermutator'] = primary_samples['SAMPLE_ID'].isin(hypermutators).astype(int)
print(f"Filtered median: {sample_counts.median():.0f}, threshold > {threshold:.0f}")
print(f"Hypermutators: {len(hypermutators)} samples")


# HPV STATUS
print("\n HPV Status")
hpv_cols = [c for c in patient_data.columns if 'HPV' in c.upper()]
print("HPV Columns:", hpv_cols)

if 'HPV_STATUS_P16' in patient_data.columns:
        patient_data['HPV_clean'] = patient_data['HPV_STATUS_P16']
else:
        patient_data['HPV_clean'] = pd.NA

print(patient_data['HPV_clean'].value_counts(dropna=False))


# FINAL MERGES
print("\n Merging...")
patient_cols =[
        'PATIENT_ID', 'AGE', 'SEX', 'PRIMARY_SITE_PATIENT',
    'AJCC_PATHOLOGIC_TUMOR_STAGE', 'OS_STATUS', 'OS_MONTHS',
    'DFS_STATUS', 'DFS_MONTHS', 'HPV_clean'
]
clinical_final = primary_samples.merge(
        patient_data[patient_cols], on='PATIENT_ID', how='left'
)

mutations_complete = mutations_final.merge(
        clinical_final[['SAMPLE_ID', 'hypermutator', 'HPV_clean', 'PRIMARY_SITE_PATIENT']], left_on='Tumor_Sample_Barcode', right_on='SAMPLE_ID', how='left'
)

print(f"Clinical: {len(primary_samples)} -> {len(clinical_final)}")
print(f"Mutations: {len(mutations_final)} -> {len(mutations_complete)}")

unlinked = mutations_complete['SAMPLE_ID'].isna().sum()
print(f"Unlinked mutations: {unlinked} ({unlinked/len(mutations_complete)*100:.1f}%)")

unlinked_samples = mutations_complete[mutations_complete['SAMPLE_ID'].isna()]['Tumor_Sample_Barcode'].unique()
print("Unlinked samples:", unlinked_samples)

# SAVE + SUMMARY
print("\n SUMMARY & SAVE")
print("_^" * 50)
print(f". Mutations:     {len(mutations):,} -> {len(mutations_final):,} coding ({100*(1-len(mutations_final)/len(mutations)):.1f}% filtered)")
print(f". Samples:       {len(sample_data)} -> {len(primary_samples)} primary")
print(f". Hypermutators: {len(hypermutators)} flagged (> {threshold:.0f})")
print(f". HPV evaluable: {clinical_final['HPV_clean'].notna().sum()}")
print(f". HPV+:          {len(mutations_complete[mutations_complete['HPV_clean']=='Positive']):,}")
print(f". Final shape:   {mutations_complete.shape}")


# SAVE
mutations_complete.to_csv('../data/clean_mutations.tsv', sep='\t', index=False)
clinical_final.to_csv('../data/clean_clinical.tsv', sep='\t', index=False)
print("\n SAVED: clean_mutations.tsv, clean_clinical.tsv")


print("COMPLETE! Ready for HPV analysis.")