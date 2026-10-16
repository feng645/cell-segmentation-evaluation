import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import json

def _row_mean(X):
    """mean over rows -> 1D numpy array; works for sparse/dense"""
    m = X.mean(axis=0)
    # m can be matrix (1, n) if sparse; flatten to 1D
    return np.asarray(m).ravel()

def read_and_aggregate_h5ad(file_path, ct_map, sample_ratio=1.0):
    """
    Read the h5ad, filter and standardize cell_type by the ct_map keys,
    then average expression per cell_type; returns a DataFrame (rows = cell_type, columns = genes).
    """
    print(f"Reading {file_path}...")
    if not file_path.endswith(".h5ad"):
        raise ValueError("Only .h5ad is supported")

    adata = sc.read_h5ad(file_path)

    # Locate the cell-type column
    if "ann_finest_level" in adata.obs.columns:
        cell_type_col = "ann_finest_level"
    else:
        raise ValueError("Cell type column not found in .h5ad")

    # Case-insensitive substring match on each key; matched names are standardized to the key
    ct_keys = list(ct_map.keys())
    # Source cell-type column as lower-case strings
    src = adata.obs[cell_type_col].astype(str).str.lower()
    # Empty Series for the matches
    matched = pd.Series(index=adata.obs.index, dtype=object)

    # Match key by key (e.g. "T", "B", "NK")
    for k in ct_keys:
        mask = src.str.contains(k.lower(), na=False) & matched.isna()
        matched.loc[mask] = k
    keep = matched.notna()
    adata = adata[keep]
    if len(adata) == 0:
        print("No cells left after ct_map filtering.")
        return pd.DataFrame()
    adata.obs["cell_type_std"] = matched.loc[adata.obs.index].values

    # Optional subsampling (before aggregation)
    if sample_ratio < 1.0 and len(adata) > 0:
        n = max(1, int(len(adata) * sample_ratio))
        idx = np.random.choice(len(adata), n, replace=False)
        adata = adata[idx]

    # Average expression per cell type
    # gene_names = adata.var_names.to_list()
    gene_names = adata.var["feature_name"].to_list()
    rows = []
    row_names = []
    for ct in sorted(adata.obs["cell_type_std"].unique()):
        ad_sub = adata[adata.obs["cell_type_std"] == ct]
        mean_vec = _row_mean(ad_sub.X)
        rows.append(mean_vec)
        row_names.append(ct)

    df = pd.DataFrame(rows, index=row_names, columns=gene_names)
    # Column names are feature names and may be duplicated
    # Average duplicated columns with groupby(axis=1)
    df = df.groupby(df.columns, axis=1).mean()

    # Add helper columns: ct_idx, cell_type, atlas
    df.insert(len(df.columns), "ct_idx", [ct_map[ct] for ct in df.index])
    df.insert(len(df.columns), "cell_type", df.index)
    df.insert(len(df.columns), "atlas", "HLCA")

    # Stable order: sort by ct_idx; columns = genes + ct_idx + cell_type + atlas
    gene_cols = [c for c in df.columns if c not in ["ct_idx", "cell_type", "atlas"]]
    df = df.loc[sorted(df.index, key=lambda x: ct_map[x]), gene_cols + ["ct_idx", "cell_type", "atlas"]]

    print(f"Aggregated: {df.shape[0]} cell types, {len(gene_cols)} genes")
    print(f"Unique cell types:\n{df[['cell_type', 'ct_idx', 'atlas']].drop_duplicates()}")
    return df

if __name__ == "__main__":
    ap = argparse.ArgumentParser("h5ad aggregator: filter by ct_map and average per cell type")
    ap.add_argument("--mode", choices=["create", "info"], default="create")
    ap.add_argument("--input_file", type=str, default="HCLA.h5ad")
    ap.add_argument("--output_file", type=str, default="lung_data_avg.csv")
    ap.add_argument("--csv_file", type=str, default="lung_data_avg.csv")
    ap.add_argument("--ct_map", type=str,
                    help='JSON dict or file path, e.g. {"T":0,"B":1,"NK":2}')
    ap.add_argument("--sample_ratio", type=float, default=1.0)
    args = ap.parse_args()

    # Parse ct_map
    ct_map = None
    if args.ct_map:
        try:
            if os.path.exists(args.ct_map):
                with open(args.ct_map, "r") as f:
                    ct_map = json.load(f)
            else:
                ct_map = json.loads(args.ct_map)
            ct_map = {str(k): int(v) for k, v in ct_map.items()}
            print(f"ct_map: {ct_map}")
        except Exception as e:
            print(f"Failed to parse --ct_map: {e}")
            ct_map = None

    if args.mode == "create":
        df = read_and_aggregate_h5ad(args.input_file, ct_map, args.sample_ratio)
        if df.empty:
            raise SystemExit(0)
        df.to_csv(args.output_file, index=True)
        print(f"Saved: {args.output_file}  shape={df.shape}")

    elif args.mode == "info":
        if not os.path.exists(args.csv_file):
            raise FileNotFoundError(f"{args.csv_file} not found.")
        df = pd.read_csv(args.csv_file, index_col=0)
        gene_cols = [c for c in df.columns if c not in ["ct_idx", "cell_type", "atlas"]]
        print("="*50)
        print(f"Dataset: {args.csv_file}")
        print("="*50)
        print(f"Pseudo-cells (cell types): {len(df)}")
        print(f"Genes: {len(gene_cols)}")
        print("\nTypes:")
        g = df[["ct_idx", "cell_type"]].drop_duplicates().sort_values("ct_idx")
        for _, r in g.iterrows():
            print(f"  [{int(r['ct_idx']):2d}] {r['cell_type']}")
