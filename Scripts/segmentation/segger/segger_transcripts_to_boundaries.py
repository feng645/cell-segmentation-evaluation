#!/usr/bin/env python3
"""
Convert Segger transcript assignments into cell boundary files.

This script is intentionally standalone so it can be inserted into the Segger
workflow without modifying the evaluation code. It reads a
`segger_transcripts.parquet` file, groups transcripts by Segger cell ID, builds
polygon boundaries, and writes outputs that the downstream R evaluation code can
consume directly.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Iterable

import pandas as pd

try:
    import geopandas as gpd
except ImportError as exc:  # pragma: no cover - workflow requires geopandas
    raise SystemExit(
        "geopandas is required for this script. Run it in the Segger conda environment."
    ) from exc

try:
    from shapely import make_valid
except ImportError:  # shapely < 2
    make_valid = None

from shapely.geometry import mapping

from segger.prediction.boundary import generate_boundaries


UNASSIGNED_VALUES = {"", "-1", "UNASSIGNED", "None", "none", "nan", "NA", "NaN"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Segger transcript parquet output into cell boundary files."
    )
    parser.add_argument(
        "--transcripts-parquet",
        required=True,
        help="Path to segger_transcripts.parquet.",
    )
    parser.add_argument(
        "--geojson-output",
        default=None,
        help="Output GeoJSON(.gz) path. Default: <input-dir>/segger_cell_boundaries.geojson.gz",
    )
    parser.add_argument(
        "--parquet-output",
        default=None,
        help="Optional GeoParquet output path.",
    )
    parser.add_argument(
        "--cell-id-col",
        default="segger_cell_id",
        help="Column containing Segger cell IDs.",
    )
    parser.add_argument(
        "--x-col",
        default=None,
        help="X coordinate column. Auto-detected if omitted.",
    )
    parser.add_argument(
        "--y-col",
        default=None,
        help="Y coordinate column. Auto-detected if omitted.",
    )
    parser.add_argument(
        "--min-points",
        type=int,
        default=3,
        help="Minimum transcripts per cell needed to attempt boundary generation.",
    )
    parser.add_argument(
        "--simplify-tolerance",
        type=float,
        default=0.0,
        help="Optional polygon simplification tolerance in coordinate units.",
    )
    return parser.parse_args()


def detect_coord_cols(df: pd.DataFrame, x_col: str | None, y_col: str | None) -> tuple[str, str]:
    if x_col and y_col:
        missing = [col for col in (x_col, y_col) if col not in df.columns]
        if missing:
            raise KeyError(f"Missing requested coordinate columns: {', '.join(missing)}")
        return x_col, y_col

    candidates = [
        ("x_location", "y_location"),
        ("global_x", "global_y"),
        ("x", "y"),
    ]
    for cand_x, cand_y in candidates:
        if cand_x in df.columns and cand_y in df.columns:
            return cand_x, cand_y

    raise KeyError(
        "Could not auto-detect coordinate columns. Tried: "
        "x_location/y_location, global_x/global_y, x/y."
    )


def normalize_cell_ids(series: pd.Series) -> pd.Series:
    out = series.astype("string")
    out = out.fillna("UNASSIGNED")
    out = out.str.strip()
    out = out.mask(out.isin(UNASSIGNED_VALUES), "UNASSIGNED")
    return out


def write_geojson(features: Iterable[dict], output_path: Path) -> None:
    payload = {"type": "FeatureCollection", "features": list(features)}
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.suffix == ".gz":
        with gzip.open(output_path, "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
    else:
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))


def gdf_to_features(gdf: gpd.GeoDataFrame) -> list[dict]:
    features: list[dict] = []
    for _, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        if make_valid is not None:
            geom = make_valid(geom)
        if geom is None or geom.is_empty:
            continue
        features.append(
            {
                "type": "Feature",
                "id": str(row["cell_id"]),
                "properties": {
                    "cell_id": str(row["cell_id"]),
                    "length": int(row["length"]),
                },
                "geometry": mapping(geom),
            }
        )
    return features


def main() -> None:
    args = parse_args()
    transcripts_path = Path(args.transcripts_parquet).expanduser().resolve()
    if not transcripts_path.exists():
        raise FileNotFoundError(f"Transcript parquet not found: {transcripts_path}")

    geojson_output = (
        Path(args.geojson_output).expanduser().resolve()
        if args.geojson_output
        else transcripts_path.with_name("segger_cell_boundaries.geojson.gz")
    )
    parquet_output = (
        Path(args.parquet_output).expanduser().resolve()
        if args.parquet_output
        else None
    )

    print(f"[load] Reading transcripts: {transcripts_path}")
    df = pd.read_parquet(transcripts_path)
    if args.cell_id_col not in df.columns:
        raise KeyError(
            f"Missing cell ID column '{args.cell_id_col}'. Available: {list(df.columns)}"
        )

    x_col, y_col = detect_coord_cols(df, args.x_col, args.y_col)
    print(f"[schema] Using cell id column: {args.cell_id_col}")
    print(f"[schema] Using coordinates: {x_col}, {y_col}")

    work = df[[args.cell_id_col, x_col, y_col]].copy()
    work[args.cell_id_col] = normalize_cell_ids(work[args.cell_id_col])
    work = work[work[args.cell_id_col] != "UNASSIGNED"]
    work = work.dropna(subset=[x_col, y_col])

    if work.empty:
        raise ValueError("No assigned transcripts remain after filtering.")

    cell_sizes = work.groupby(args.cell_id_col).size()
    keep_ids = cell_sizes[cell_sizes >= args.min_points].index
    work = work[work[args.cell_id_col].isin(keep_ids)].copy()

    print(f"[filter] Kept {len(keep_ids)} cells with >= {args.min_points} transcripts")
    print(f"[filter] Retained {len(work):,} assigned transcripts")

    boundaries = generate_boundaries(
        work,
        x=x_col,
        y=y_col,
        cell_id=args.cell_id_col,
    )
    boundaries = boundaries.dropna(subset=["geometry"]).copy()
    boundaries = boundaries[~boundaries.geometry.is_empty].copy()

    if args.simplify_tolerance > 0:
        boundaries["geometry"] = boundaries.geometry.simplify(
            tolerance=args.simplify_tolerance,
            preserve_topology=True,
        )
        boundaries = boundaries.dropna(subset=["geometry"]).copy()
        boundaries = boundaries[~boundaries.geometry.is_empty].copy()

    boundaries["cell_id"] = boundaries["cell_id"].astype(str)
    boundaries["length"] = boundaries["length"].astype(int)

    features = gdf_to_features(boundaries)
    if not features:
        raise ValueError("Boundary generation produced no valid geometries.")

    print(f"[write] Writing GeoJSON: {geojson_output}")
    write_geojson(features, geojson_output)

    if parquet_output is not None:
        parquet_output.parent.mkdir(parents=True, exist_ok=True)
        print(f"[write] Writing GeoParquet: {parquet_output}")
        boundaries.to_parquet(parquet_output)

    print(
        f"[done] Generated {len(features)} cell boundaries from "
        f"{len(work):,} transcripts."
    )


if __name__ == "__main__":
    main()
