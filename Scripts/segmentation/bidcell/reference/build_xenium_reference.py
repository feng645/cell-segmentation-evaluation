#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


HCLA_CT_MAP = {
    "AT2": 0,
    "AT1": 1,
    "Smooth muscle": 2,
    "CD8 T cells": 3,
    "CD4 T cells": 4,
    "T cells proliferating": 5,
    "NK cells": 6,
    "B cells": 7,
    "Plasma cells": 8,
    "Classical monocytes": 9,
    "Alveolar macrophages": 10,
    "Monocyte-derived Mph": 11,
    "Fibroblasts": 12,
    "Endothelial cells": 13,
}

NSCLC_CT_MAP = {
    "AT2": 0,
    "AT1": 1,
    "Smooth muscle": 2,
    "CD8 T cells": 3,
    "CD4 T cells": 4,
    "T cells proliferating": 5,
    "NK cells": 6,
    "B cells": 7,
    "Plasma cells": 8,
    "Classical monocytes": 9,
    "Alveolar macrophages": 10,
    "Monocyte-derived Mph": 11,
    "Fibroblasts": 12,
    "Endothelial": 13,
    "Malignant": 14,
}

CELL_TYPE_SORT_ORDER = dict(HCLA_CT_MAP)
CELL_TYPE_SORT_ORDER["Endothelial"] = 13
CELL_TYPE_SORT_ORDER["Malignant"] = 14

CONTROL_PREFIXES = [
    "NegControlProbe_",
    "NegControlCodeword_",
    "antisense_",
    "BLANK_",
    "Blank-",
    "NegPrb",
]


def row_mean(x):
    mean_value = x.mean(axis=0)
    return np.asarray(mean_value).ravel()


def is_control_feature(feature_name, control_prefixes):
    feature_name = str(feature_name)
    return any(feature_name.startswith(prefix) for prefix in control_prefixes)


def read_xenium_panel(transcripts_fp, gene_col, control_prefixes):
    panel = pd.read_csv(transcripts_fp, usecols=[gene_col])[gene_col].dropna().astype(str)
    panel_unique = sorted(set(panel.tolist()))
    panel_clean = [gene for gene in panel_unique if not is_control_feature(gene, control_prefixes)]
    return panel_unique, panel_clean


def standardize_hcla_types(obs_series):
    matched = pd.Series(index=obs_series.index, dtype=object)
    src = obs_series.astype(str).str.lower()

    for cell_type in HCLA_CT_MAP.keys():
        mask = src.str.contains(cell_type.lower(), na=False) & matched.isna()
        matched.loc[mask] = cell_type

    return matched


def build_hcla_reference(hcla_fp):
    adata = sc.read_h5ad(hcla_fp)

    if "ann_finest_level" not in adata.obs.columns:
        raise ValueError("Expected 'ann_finest_level' in HCLA .h5ad obs.")

    matched = standardize_hcla_types(adata.obs["ann_finest_level"])
    keep = matched.notna()
    adata = adata[keep].copy()

    if adata.n_obs == 0:
        raise ValueError("No HCLA cells matched the configured target cell types.")

    adata.obs["cell_type_std"] = matched.loc[adata.obs.index].values

    if "feature_name" in adata.var.columns:
        gene_names = adata.var["feature_name"].astype(str).tolist()
    else:
        gene_names = adata.var_names.astype(str).tolist()

    rows = []
    row_names = []
    for cell_type in sorted(adata.obs["cell_type_std"].unique(), key=lambda x: HCLA_CT_MAP[x]):
        ad_sub = adata[adata.obs["cell_type_std"] == cell_type]
        rows.append(row_mean(ad_sub.X))
        row_names.append(cell_type)

    df = pd.DataFrame(rows, index=row_names, columns=gene_names)
    df = df.T.groupby(level=0).mean().T
    df.insert(len(df.columns), "ct_idx", [HCLA_CT_MAP[cell_type] for cell_type in df.index])
    df.insert(len(df.columns), "cell_type", df.index)
    df.insert(len(df.columns), "atlas", "HLCA")
    return df.reset_index(drop=True)


def detect_sep(txt_fp):
    with open(txt_fp, "r", encoding="utf-8") as handle:
        preview = "".join([next(handle) for _ in range(3)])
    return "\t" if "\t" in preview else r"\s+"


def build_nsclc_reference(txt_fp):
    sep = detect_sep(txt_fp)
    df = pd.read_csv(txt_fp, sep=sep, engine="python", header=0, index_col=0)
    df = df.groupby(level=0).mean().T

    ct_keys = list(NSCLC_CT_MAP.keys())
    df = df[df.index.isin(ct_keys)].copy()

    if df.empty:
        print(f"[warn] no supported NSCLC columns found in {txt_fp}")
        return pd.DataFrame()

    df.insert(len(df.columns), "ct_idx", [NSCLC_CT_MAP[cell_type] for cell_type in df.index])
    df.insert(len(df.columns), "cell_type", df.index)
    df.insert(len(df.columns), "atlas", Path(txt_fp).stem.replace("_expression_Celltype_majorlineage", ""))
    return df.reset_index(drop=True)


def merge_references(reference_frames):
    frames = [frame.copy() for frame in reference_frames if not frame.empty]
    if not frames:
        raise ValueError("No reference tables were generated.")

    metadata_cols = ["ct_idx", "cell_type", "atlas"]
    common_genes = set([col for col in frames[0].columns if col not in metadata_cols])
    for frame in frames[1:]:
        genes = set([col for col in frame.columns if col not in metadata_cols])
        common_genes &= genes

    common_genes = sorted(common_genes)
    if not common_genes:
        raise ValueError("No common genes found across the reference sources.")

    merged = pd.concat([frame[common_genes + metadata_cols] for frame in frames], axis=0, ignore_index=True)
    return merged, common_genes


def subset_and_scale_reference(reference_df, panel_genes):
    metadata_cols = ["ct_idx", "cell_type", "atlas"]
    gene_cols = [gene for gene in panel_genes if gene in reference_df.columns]

    if not gene_cols:
        raise ValueError("No Xenium panel genes overlap the merged reference.")

    final_df = reference_df[gene_cols + metadata_cols].copy()

    values = final_df[gene_cols].astype(float)
    row_min = values.min(axis=1)
    row_max = values.max(axis=1)
    row_range = (row_max - row_min).replace(0, 1)
    final_df.loc[:, gene_cols] = values.sub(row_min, axis=0).div(row_range, axis=0)

    return final_df, gene_cols


def derive_markers(reference_df, max_overlaps_pos, max_overlaps_neg):
    n_genes = reference_df.shape[1] - 3
    ref_expr = reference_df.iloc[:, :n_genes].to_numpy(dtype=float)
    gene_names = reference_df.columns[:n_genes]

    cell_type_table = reference_df[["ct_idx", "cell_type"]].drop_duplicates().sort_values(["ct_idx", "cell_type"])
    cell_types = cell_type_table["cell_type"].tolist()
    n_cell_types = len(cell_types)

    ct_idx = reference_df["ct_idx"].to_numpy()

    pct_10 = np.percentile(ref_expr, 10, axis=1, keepdims=True)
    low_expr_true = (ref_expr <= pct_10).astype(int)
    low_expr_true_agg = np.zeros((n_cell_types, n_genes))
    for cell_type in cell_types:
        idx = CELL_TYPE_SORT_ORDER[cell_type]
        rows = np.where(ct_idx == idx)[0]
        low_expr_true_agg[cell_types.index(cell_type), :] = np.prod(low_expr_true[rows], axis=0)
    overlaps_neg = np.sum(low_expr_true_agg, axis=0)
    low_expr_true_agg[:, overlaps_neg > max_overlaps_neg] = 0
    df_neg = pd.DataFrame(low_expr_true_agg, index=cell_types, columns=gene_names)

    pct_90 = np.percentile(ref_expr, 90, axis=1, keepdims=True)
    high_expr_true = (ref_expr >= pct_90).astype(int)
    high_expr_true_agg = np.zeros((n_cell_types, n_genes))
    for cell_type in cell_types:
        idx = CELL_TYPE_SORT_ORDER[cell_type]
        rows = np.where(ct_idx == idx)[0]
        high_expr_true_agg[cell_types.index(cell_type), :] = np.prod(high_expr_true[rows], axis=0)
    overlaps_pos = np.sum(high_expr_true_agg, axis=0)
    high_expr_true_agg[:, overlaps_pos > max_overlaps_pos] = 0
    df_pos = pd.DataFrame(high_expr_true_agg, index=cell_types, columns=gene_names)

    return df_pos, df_neg


def write_lines(output_fp, values):
    with open(output_fp, "w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Build Xenium-compatible BIDCell lung reference and marker files.")
    parser.add_argument(
        "--transcripts",
        default="/scratch.global/feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts.csv.gz",
        help="Xenium transcripts.csv.gz path.",
    )
    parser.add_argument(
        "--hcla",
        default="/users/9/feng0421/BIDCell/data/sc_references/HCLA.h5ad",
        help="HCLA .h5ad path.",
    )
    parser.add_argument(
        "--nsclc",
        nargs="+",
        default=[
            "/users/9/feng0421/BIDCell/data/sc_references/NSCLC_GSE117570_expression_Celltype_majorlineage.txt",
            "/users/9/feng0421/BIDCell/data/sc_references/NSCLC_GSE143423_expression_Celltype_majorlineage.txt",
            "/users/9/feng0421/BIDCell/data/sc_references/NSCLC_GSE148071_expression_Celltype_majorlineage.txt",
        ],
        help="One or more NSCLC atlas txt files.",
    )
    parser.add_argument(
        "--output-dir",
        default="/users/9/feng0421/BIDCell/data/sc_references",
        help="Directory for output reference files.",
    )
    parser.add_argument(
        "--gene-col",
        default="feature_name",
        help="Gene column in Xenium transcripts.",
    )
    parser.add_argument(
        "--max-overlaps-pos",
        type=int,
        default=4,
        help="Maximum positive-marker overlap across cell types.",
    )
    parser.add_argument(
        "--max-overlaps-neg",
        type=int,
        default=15,
        help="Maximum negative-marker overlap across cell types.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    panel_all, panel_clean = read_xenium_panel(
        transcripts_fp=args.transcripts,
        gene_col=args.gene_col,
        control_prefixes=CONTROL_PREFIXES,
    )
    print(f"[info] Xenium panel genes: total={len(panel_all)} clean={len(panel_clean)}")

    hcla_df = build_hcla_reference(args.hcla)
    print(f"[info] HCLA reference rows={len(hcla_df)}")

    nsclc_frames = []
    for nsclc_fp in args.nsclc:
        frame = build_nsclc_reference(nsclc_fp)
        if not frame.empty:
            print(f"[info] NSCLC reference rows={len(frame)} from {nsclc_fp}")
            nsclc_frames.append(frame)

    merged_df, merged_genes = merge_references([hcla_df] + nsclc_frames)
    print(f"[info] merged reference rows={len(merged_df)} common_genes={len(merged_genes)}")

    final_ref_df, final_gene_cols = subset_and_scale_reference(merged_df, panel_clean)
    print(f"[info] final Xenium reference rows={len(final_ref_df)} genes={len(final_gene_cols)}")

    pos_df, neg_df = derive_markers(
        reference_df=final_ref_df,
        max_overlaps_pos=args.max_overlaps_pos,
        max_overlaps_neg=args.max_overlaps_neg,
    )

    genes_in_reference = set(final_gene_cols)
    targets_to_filter = sorted([gene for gene in panel_all if gene not in genes_in_reference])

    raw_ref_fp = output_dir / "lung_ref_xenium_raw.csv"
    final_ref_fp = output_dir / "lung_ref_xenium.csv"
    pos_fp = output_dir / "xenium_lung_markers_pos.csv"
    neg_fp = output_dir / "xenium_lung_markers_neg.csv"
    panel_fp = output_dir / "xenium_gene_panel.txt"
    common_fp = output_dir / "xenium_common_genes.txt"
    filter_fp = output_dir / "xenium_target_to_filter.txt"

    merged_df.to_csv(raw_ref_fp, index=True)
    final_ref_df.to_csv(final_ref_fp, index=True)
    pos_df.to_csv(pos_fp)
    neg_df.to_csv(neg_fp)
    write_lines(panel_fp, panel_clean)
    write_lines(common_fp, merged_genes)
    write_lines(filter_fp, targets_to_filter)

    print("[done] wrote:")
    for fp in [raw_ref_fp, final_ref_fp, pos_fp, neg_fp, panel_fp, common_fp, filter_fp]:
        print(f"  {fp}")


if __name__ == "__main__":
    main()
