# so far, you should already get a lung_ref.csv and run get_gene_list.py to get gene list txt file

import pandas as pd

# Input and output paths
ref_csv = "/scratch.global/feng0421/data/ref/lung_ref.csv"
gene_txt = "/scratch.global/feng0421/data/MERFISH/dataset/Merfish_gene.txt"
out_csv = "/scratch.global/feng0421/data/ref/lung_ref_merfish.csv"

# Read gene names to keep
with open(gene_txt) as f:
    keep_genes = set(line.strip() for line in f if line.strip())

# Load reference file
df = pd.read_csv(ref_csv)

# Define the slicing range: from 2nd column to the 4th-from-last
start = 1                       # second column (0-based index)
end = len(df.columns) - 3       # exclude last 3 columns
gene_cols = df.columns[start:end]

# Select only those gene columns that are in the keep list
gene_cols_keep = [c for c in gene_cols if c in keep_genes]

# Build the final DataFrame without changing column order
# Start: first column(s), then filtered genes, then last 3 columns

final_columns = [df.columns[0]] + gene_cols_keep + df.columns[end:].tolist()

df_filtered = df[final_columns]

# Scaling gene columns by row (MinMax scaling to [0, 1])
print("Scaling gene columns by row (Min-Max 0-1), skipping header row...")

# Gene columns: from the 2nd to the 4th-from-last column
gene_start_idx = 1  # starts at the 2nd column (AATK)
gene_end_idx = len(df_filtered.columns) - 3  # exclude the last 3 columns

# Check the number of rows
if len(df_filtered) > 1:  # make sure there are rows to process
    # Scale rows from the 2nd onward (skip the header)
    gene_data = df_filtered.iloc[1:, gene_start_idx:gene_end_idx]
    
    # Make sure the data are numeric
    gene_data = gene_data.astype(float)
    
    # Row-wise min-max scaling
    row_min = gene_data.min(axis=1)
    row_max = gene_data.max(axis=1)
    row_range = (row_max - row_min).replace(0, 1)  # avoid division by zero
    
    gene_data_scaled = gene_data.sub(row_min, axis=0).div(row_range, axis=0)
    
    # Put the scaled values back into the DataFrame
    df_filtered.iloc[1:, gene_start_idx:gene_end_idx] = gene_data_scaled
    
    print(f"Scaled {len(gene_data)} rows and {len(gene_data.columns)} gene columns")
else:
    print("No data rows to scale (only header row present)")

# Check for missing genes
missing_genes = list(keep_genes - set(gene_cols))
print(f"gene in gene.txt but not in ref file: {missing_genes}")
print(f"Total genes missing: {len(missing_genes)}")

# Save
df_filtered = df_filtered.reset_index(drop=True)
df_filtered = df_filtered.drop(df_filtered.columns[0], axis=1)
df_filtered.to_csv(out_csv, index=True)
print(f"Filtered reference saved: {out_csv}")
print(f"Kept {len(gene_cols_keep)} genes")
