#!/usr/bin/env python3
"""

Load cBioPortal mutation and clinical data into DataFrames.


"""

import pandas as pd
import os


# File paths
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

# Print info
for name, df in [('MUTATIONS', mutations), ('PATIENT',
                                            patient_data),('SAMPLE',sample_data)]:
    print(f"\n=== {name} ===")
    print(f"Shape: {df.shape}")
    print("Columns:", list(df.columns[:8]), "...")
    print(df.head(2))


# HPV status investigation
for col in ['HPV_STATUS_P16', 'HPV_STATUS_ISH']:
    if col in patient_data.columns:
        print(f"\n=== {col} ===")
        print(f"Non-null: {patient_data[col].notna().sum()} / {len(patient_data)}")
        print(patient_data[col].value_counts(dropna=False))


# NOTE: clinical data uses "[Not Available]" and "[Not Evaluated]" strings
# instead of NaN for missing values. Must handle during cleaning step.
# HPV decision: Use HPV_STATUS_P16 (n=115) over HPV_STATUS_ISH (n=89) for
# better statistical power in stratified analysis.
