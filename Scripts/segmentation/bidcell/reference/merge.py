#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def main(file1, file2, output, common_txt):
    # Read the two CSV files, using the first column as index
    df1 = pd.read_csv(file1, index_col=0)
    df2 = pd.read_csv(file2, index_col=0)

    # Find common columns (genes)
    # Keep the order from df1 to avoid shuffling
    common_cols = [col for col in df1.columns if col in df2.columns]

    if not common_cols:
        raise ValueError("No common columns found between the two files.")

    # Save the list of common columns into a txt file
    with open(common_txt, "w") as f:
        for col in common_cols:
            f.write(col + "\n")
    print(f"Common columns saved to: {common_txt} (count={len(common_cols)})")

    # Subset both DataFrames to only the common columns
    df1_sub = df1[common_cols]
    df2_sub = df2[common_cols]

    # Concatenate along rows, keeping original order
    merged = pd.concat([df1_sub, df2_sub], axis=0)

    # Reset index to 0,1,2...
    merged = merged.reset_index(drop=True)
    merged = merged.drop(merged.columns[0], axis=1)
    # Save merged CSV
    merged.to_csv(output, index=True)
    print(f"Merged file saved: {output}")
    print(f"Shape: {merged.shape} (rows: {merged.shape[0]}, columns: {merged.shape[1]})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge two CSV files by common columns")
    parser.add_argument("--file1", required=True, help="Path to the first CSV file")
    parser.add_argument("--file2", required=True, help="Path to the second CSV file")
    parser.add_argument("--output", default="merged.csv", help="Path to save the merged CSV file")
    parser.add_argument("--common_txt", default="common_columns.txt", help="Path to save the common columns list")
    args = parser.parse_args()

    main(args.file1, args.file2, args.output, args.common_txt)
# Example usage:
# cd /users/9/feng0421/Code/BIDCell/sc_ref
# python merge.py --file1 /users/9/feng0421/Code/BIDCell/sc_ref/lung_data.csv --file2 /scratch.global/feng0421/data/ref/NSCLC_GSE148071_expression_Celltype_majorlineage.csv --output /scratch.global/feng0421/data/ref/lung_ref.csv --common_txt common_genes.txt
# python merge.py --file1 /scratch.global/feng0421/data/ref/lung_ref.csv --file2 /scratch.global/feng0421/data/ref/NSCLC_GSE117570_expression_Celltype_majorlineage.csv --output /scratch.global/feng0421/data/ref/lung_ref.csv --common_txt common_genes.txt