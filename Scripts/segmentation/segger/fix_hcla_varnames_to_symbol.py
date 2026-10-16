#!/usr/bin/env python3
"""
Create a symbol-based H5AD file for Segger.

This script rewrites `var_names` to `var[feature_name]` so Segger can match
MERFISH gene symbols instead of Ensembl IDs in `var_names`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import scanpy as sc


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Rewrite AnnData var_names to gene symbols for Segger."
    )
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        required=True,
        help="Input h5ad path (e.g. HCLA.h5ad).",
    )
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        required=True,
        help="Output h5ad path with symbol-based var_names.",
    )
    parser.add_argument(
        "--symbol-column",
        type=str,
        default="feature_name",
        help="Column in adata.var used as gene symbols (default: feature_name).",
    )
    parser.add_argument(
        "--allow-overwrite",
        action="store_true",
        help="Overwrite output file if it already exists.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    in_path = args.input_h5ad
    out_path = args.output_h5ad

    if not in_path.exists():
        print(f"[ERROR] Input file not found: {in_path}", file=sys.stderr)
        return 1

    if out_path.exists() and not args.allow_overwrite:
        print(
            f"[ERROR] Output already exists: {out_path}\n"
            "Use --allow-overwrite to replace it.",
            file=sys.stderr,
        )
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Reading: {in_path}")
    adata = sc.read_h5ad(in_path)

    if args.symbol_column not in adata.var.columns:
        print(
            f"[ERROR] Column '{args.symbol_column}' not found in adata.var.\n"
            f"Available columns: {list(adata.var.columns)}",
            file=sys.stderr,
        )
        return 1

    old_preview = [str(x) for x in adata.var_names[:10]]
    adata.var_names = adata.var[args.symbol_column].astype(str)
    # Avoid AnnData write error when index.name collides with an existing
    # var column name (e.g. index.name == "feature_name").
    adata.var_names.name = "gene_symbol"
    adata.var_names_make_unique()
    new_preview = [str(x) for x in adata.var_names[:10]]

    print(f"[INFO] n_vars: {adata.n_vars}")
    print(f"[INFO] Old var_names preview: {old_preview}")
    print(f"[INFO] New var_names preview: {new_preview}")
    print(f"[INFO] Writing: {out_path}")

    adata.write_h5ad(out_path)
    print("[INFO] Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
