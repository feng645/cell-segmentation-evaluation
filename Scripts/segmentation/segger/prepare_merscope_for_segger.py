#!/usr/bin/env python3
"""
Prepare MERSCOPE raw exports into Segger-compatible parquet files.

Outputs (under --output-dir):
1) detected_transcripts.parquet
2) cellpose_nucleus_micron_space.parquet

This version derives proxy cell boundaries directly from transcript-level
cell assignments in the raw CSV. The generated boundary GeoParquet preserves
its index as EntityID so Segger can recover positive tx->bd edges even when
upstream boundary files use incompatible indexing.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from shapely.geometry import box


REQUIRED_TX_COLUMNS = [
    "global_x",
    "global_y",
    "global_z",
    "gene",
    "cell_id",
    "EntityID",
    "overlaps_nucleus",
    "transcript_id_original",
    "transcript_id",
]

REQUIRED_BD_COLUMNS = ["Geometry", "EntityID", "centroid_x", "centroid_y"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert MERSCOPE raw data into Segger parquet inputs."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("/users/9/feng0421/BIDCell/data/MerFish"),
        help="MERSCOPE data directory containing the transcript CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for Segger-ready parquet files. Default: <input-dir>/segger_ready",
    )
    parser.add_argument(
        "--transcripts-csv",
        type=Path,
        default=None,
        help="Path to transcript CSV. If omitted, auto-detects preferred file.",
    )
    parser.add_argument(
        "--boundaries-dir",
        type=Path,
        default=None,
        help="Legacy option. Ignored because boundaries are derived from transcript assignments.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=1_000_000,
        help="CSV/parquet batch size for transcript conversion and aggregation.",
    )
    parser.add_argument(
        "--z-index",
        type=int,
        default=3,
        help="Legacy option. Ignored because HDF5 boundaries are no longer read.",
    )
    parser.add_argument(
        "--max-boundary-files",
        type=int,
        default=None,
        help="Legacy option. Ignored because HDF5 boundaries are no longer read.",
    )
    parser.add_argument(
        "--min-boundary-span",
        type=float,
        default=1.0,
        help="Minimum x/y span in microns for generated proxy boundary boxes.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output parquet files.",
    )
    return parser.parse_args()


def choose_transcript_csv(input_dir: Path) -> Path:
    candidates = [
        input_dir / "detected_transcripts_seg_7.csv",
        input_dir / "HumanLungCancerPatient1_detected_transcripts.csv",
        input_dir / "HumanLungCancerPatient1_detected_transcripts_filtered.csv",
        input_dir / "roi_transcripts.csv",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        "No transcript CSV candidate found. Provide --transcripts-csv explicitly."
    )


def ensure_columns(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    if "Unnamed: 0" in work.columns:
        work = work.drop(columns=["Unnamed: 0"])

    missing_xyz = [c for c in ("global_x", "global_y", "global_z") if c not in work.columns]
    if missing_xyz:
        raise ValueError(f"Missing required coordinate columns: {missing_xyz}")
    if "gene" not in work.columns:
        raise ValueError("Missing required transcript column: gene")
    if "cell_id" not in work.columns:
        raise ValueError(
            "Missing required assignment column: cell_id. "
            "This preparation workflow derives boundaries directly from transcript assignments."
        )

    work["cell_id"] = pd.to_numeric(work["cell_id"], errors="coerce").fillna(-1).astype("int64")
    work["EntityID"] = work["cell_id"]
    work["overlaps_nucleus"] = (work["cell_id"] != -1).astype("int8")

    if "transcript_id" not in work.columns:
        raise ValueError(
            "Missing required source transcript column: transcript_id. "
            "This workflow preserves the raw transcript ID for traceability."
        )
    work["transcript_id_original"] = work["transcript_id"].astype("string")

    # Rebuild transcript_id as a contiguous numeric key for Segger while keeping
    # the original source identifier in transcript_id_original.
    work["transcript_id"] = np.arange(len(work), dtype=np.int64)

    return work[REQUIRED_TX_COLUMNS].copy()


def convert_transcripts_csv_to_parquet(
    csv_path: Path, parquet_path: Path, chunksize: int, overwrite: bool
) -> Tuple[int, int]:
    if parquet_path.exists() and not overwrite:
        print(f"[transcripts] Exists, skip: {parquet_path}")
        meta = pq.read_metadata(parquet_path)
        return meta.num_rows, meta.num_columns

    if parquet_path.exists():
        parquet_path.unlink()

    print(f"[transcripts] Converting: {csv_path}")
    print(f"[transcripts] Writing to: {parquet_path}")

    writer: Optional[pq.ParquetWriter] = None
    total_rows = 0
    id_offset = 0

    for idx, chunk in enumerate(pd.read_csv(csv_path, chunksize=chunksize)):
        cleaned = ensure_columns(chunk)
        cleaned["transcript_id"] = np.arange(
            id_offset, id_offset + len(cleaned), dtype=np.int64
        )

        id_offset += len(cleaned)
        total_rows += len(cleaned)

        table = pa.Table.from_pandas(cleaned, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(parquet_path, table.schema, compression="snappy")
        writer.write_table(table)

        if (idx + 1) % 10 == 0:
            print(f"[transcripts] Processed chunks: {idx + 1}, rows: {total_rows:,}")

    if writer is None:
        raise RuntimeError("No transcript rows were read from CSV.")
    writer.close()

    print(f"[transcripts] Done. Total rows: {total_rows:,}")
    return total_rows, len(REQUIRED_TX_COLUMNS)


def summarize_assigned_cells(
    transcripts_path: Path, chunksize: int
) -> Tuple[Dict[int, list], int, int]:
    print("[boundaries] Aggregating transcript assignments into proxy cell boxes...")

    pf = pq.ParquetFile(transcripts_path)
    stats: Dict[int, list] = {}
    assigned_rows = 0
    total_rows = 0

    for idx, batch in enumerate(
        pf.iter_batches(batch_size=chunksize, columns=["global_x", "global_y", "cell_id"]),
        start=1,
    ):
        chunk = batch.to_pandas()
        if chunk.empty:
            continue

        chunk["cell_id"] = pd.to_numeric(chunk["cell_id"], errors="coerce").fillna(-1).astype("int64")
        assigned = chunk.loc[chunk["cell_id"] != -1, ["cell_id", "global_x", "global_y"]].copy()
        total_rows += len(chunk)
        assigned_rows += len(assigned)

        if assigned.empty:
            continue

        grouped = assigned.groupby("cell_id", sort=False).agg(
            min_x=("global_x", "min"),
            max_x=("global_x", "max"),
            min_y=("global_y", "min"),
            max_y=("global_y", "max"),
            sum_x=("global_x", "sum"),
            sum_y=("global_y", "sum"),
            n=("cell_id", "size"),
        )

        for entity_id, row in grouped.iterrows():
            vals = [
                float(row["min_x"]),
                float(row["max_x"]),
                float(row["min_y"]),
                float(row["max_y"]),
                float(row["sum_x"]),
                float(row["sum_y"]),
                int(row["n"]),
            ]
            prev = stats.get(int(entity_id))
            if prev is None:
                stats[int(entity_id)] = vals
            else:
                prev[0] = min(prev[0], vals[0])
                prev[1] = max(prev[1], vals[1])
                prev[2] = min(prev[2], vals[2])
                prev[3] = max(prev[3], vals[3])
                prev[4] += vals[4]
                prev[5] += vals[5]
                prev[6] += vals[6]

        if idx % 10 == 0:
            print(
                f"[boundaries] Processed batches: {idx}, rows: {total_rows:,}, assigned: {assigned_rows:,}, cells: {len(stats):,}"
            )

    if assigned_rows == 0 or not stats:
        raise RuntimeError(
            "[boundaries] No assigned transcripts were found in cell_id. "
            "Cannot derive proxy boundaries."
        )

    unassigned_rows = total_rows - assigned_rows
    print(
        f"[boundaries] Aggregation done. total_rows={total_rows:,}, assigned_rows={assigned_rows:,}, "
        f"unassigned_rows={unassigned_rows:,}, cells={len(stats):,}"
    )
    return stats, assigned_rows, total_rows


def build_assignment_boundaries(
    cell_stats: Dict[int, list],
    out_parquet: Path,
    min_boundary_span: float,
    overwrite: bool,
) -> Tuple[int, int]:
    if out_parquet.exists() and not overwrite:
        print(f"[boundaries] Exists, skip: {out_parquet}")
        gdf = gpd.read_parquet(out_parquet)
        return len(gdf), len(gdf.columns)

    if out_parquet.exists():
        out_parquet.unlink()

    half_span = max(float(min_boundary_span), 1e-6) / 2.0
    records = []

    for entity_id in sorted(cell_stats):
        min_x, max_x, min_y, max_y, sum_x, sum_y, n_tx = cell_stats[entity_id]
        centroid_x = sum_x / n_tx
        centroid_y = sum_y / n_tx

        if (max_x - min_x) < (2.0 * half_span):
            min_x = centroid_x - half_span
            max_x = centroid_x + half_span
        if (max_y - min_y) < (2.0 * half_span):
            min_y = centroid_y - half_span
            max_y = centroid_y + half_span

        records.append(
            {
                "Geometry": box(min_x, min_y, max_x, max_y),
                "EntityID": int(entity_id),
                "centroid_x": float(centroid_x),
                "centroid_y": float(centroid_y),
                "n_transcripts": int(n_tx),
            }
        )

    if not records:
        raise RuntimeError("[boundaries] No proxy boundary geometries were generated.")

    gdf = gpd.GeoDataFrame(records, geometry="Geometry")
    gdf["EntityID"] = pd.to_numeric(gdf["EntityID"], errors="raise").astype("int64")

    # Preserve EntityID as the GeoDataFrame index because Segger currently uses
    # boundary.index, not only boundary['EntityID'], when creating tx->bd edges.
    gdf.index = pd.Index(gdf["EntityID"].to_numpy(dtype=np.int64), name="boundary_index")
    gdf = gdf[["Geometry", "EntityID", "centroid_x", "centroid_y", "n_transcripts"]]
    gdf.to_parquet(out_parquet, index=True)

    print(f"[boundaries] Done. Cells: {len(gdf):,} -> {out_parquet}")
    return len(gdf), len(gdf.columns)


def validate_outputs(transcripts_path: Path, boundaries_path: Path, chunksize: int) -> None:
    print("[validate] Checking output parquet schemas and key invariants...")

    pf = pq.ParquetFile(transcripts_path)
    tx_columns = pf.schema_arrow.names
    missing_tx = [c for c in REQUIRED_TX_COLUMNS if c not in tx_columns]
    if missing_tx:
        raise ValueError(f"[validate] Missing transcript columns: {missing_tx}")

    total_rows = 0
    assigned_rows = 0
    expected_start = 0
    assigned_entity_ids = set()

    for batch in pf.iter_batches(
        batch_size=chunksize,
        columns=["cell_id", "EntityID", "transcript_id", "transcript_id_original"],
    ):
        chunk = batch.to_pandas()
        if chunk.empty:
            continue

        cell_id = pd.to_numeric(chunk["cell_id"], errors="coerce")
        entity_id = pd.to_numeric(chunk["EntityID"], errors="coerce")
        tx_id = pd.to_numeric(chunk["transcript_id"], errors="coerce")
        if cell_id.isna().any() or entity_id.isna().any() or tx_id.isna().any():
            raise ValueError("[validate] Transcript parquet contains NaN/non-numeric cell/entity/transcript_id values.")
        if chunk["transcript_id_original"].isna().any():
            raise ValueError("[validate] transcript_id_original contains missing values.")

        cell_arr = cell_id.astype("int64").to_numpy()
        entity_arr = entity_id.astype("int64").to_numpy()
        tx_arr = tx_id.astype("int64").to_numpy()

        if not np.array_equal(cell_arr, entity_arr):
            raise ValueError("[validate] cell_id and EntityID must be identical in transcripts.")

        expected = np.arange(expected_start, expected_start + len(tx_arr), dtype=np.int64)
        if not np.array_equal(tx_arr, expected):
            raise ValueError(
                "[validate] transcript_id is not the expected contiguous 0..N-1 sequence."
            )

        assigned_mask = cell_arr != -1
        assigned_rows += int(assigned_mask.sum())
        assigned_entity_ids.update(cell_arr[assigned_mask].tolist())
        total_rows += len(tx_arr)
        expected_start += len(tx_arr)

    if total_rows == 0:
        raise ValueError("[validate] Transcript parquet is empty.")
    if assigned_rows == 0:
        raise ValueError("[validate] No assigned transcripts found in cell_id/EntityID.")

    bd = gpd.read_parquet(boundaries_path)
    missing_bd = [c for c in REQUIRED_BD_COLUMNS if c not in bd.columns]
    if missing_bd:
        raise ValueError(f"[validate] Missing boundary columns: {missing_bd}")
    if bd.empty:
        raise ValueError("[validate] Boundary parquet is empty.")

    bd_entity = pd.to_numeric(bd["EntityID"], errors="coerce")
    if bd_entity.isna().any():
        raise ValueError("[validate] Boundary EntityID contains NaN/non-numeric values.")
    bd_entity = bd_entity.astype("int64")
    if bd_entity.nunique() != len(bd_entity):
        raise ValueError("[validate] Boundary EntityID values must be unique.")

    bd_index = pd.Index(bd.index)
    if len(bd_index) != len(bd_entity) or not np.array_equal(bd_index.to_numpy(dtype=np.int64), bd_entity.to_numpy(dtype=np.int64)):
        raise ValueError(
            "[validate] Boundary GeoParquet index must equal EntityID to match Segger's current reader behavior."
        )

    missing_boundaries = assigned_entity_ids.difference(set(bd_entity.tolist()))
    if missing_boundaries:
        first_missing = sorted(missing_boundaries)[:10]
        raise ValueError(
            f"[validate] Some assigned transcript EntityID values are missing from boundary parquet: {first_missing}"
        )

    print("[validate] OK")
    print(f"[validate] transcripts rows: {total_rows:,}, assigned_rows: {assigned_rows:,}")
    print(f"[validate] boundaries rows: {len(bd):,}, columns: {list(bd.columns)}")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = (args.output_dir or (input_dir / "segger_ready")).resolve()
    boundaries_dir = (args.boundaries_dir or (input_dir / "cell_boundaries")).resolve()
    transcripts_csv = args.transcripts_csv.resolve() if args.transcripts_csv else choose_transcript_csv(input_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    tx_out = output_dir / "detected_transcripts.parquet"
    bd_out = output_dir / "cellpose_nucleus_micron_space.parquet"

    print(f"[config] input_dir={input_dir}")
    print(f"[config] transcripts_csv={transcripts_csv}")
    print(f"[config] output_dir={output_dir}")
    print(f"[config] boundaries_dir={boundaries_dir} (ignored; using transcript assignments instead)")
    if args.boundaries_dir is not None or args.max_boundary_files is not None or args.z_index != 3:
        print("[config] legacy boundary options detected; they are ignored in assignment-derived mode.")

    convert_transcripts_csv_to_parquet(
        csv_path=transcripts_csv,
        parquet_path=tx_out,
        chunksize=args.chunksize,
        overwrite=args.overwrite,
    )
    cell_stats, _, _ = summarize_assigned_cells(
        transcripts_path=tx_out,
        chunksize=args.chunksize,
    )
    build_assignment_boundaries(
        cell_stats=cell_stats,
        out_parquet=bd_out,
        min_boundary_span=args.min_boundary_span,
        overwrite=args.overwrite,
    )
    validate_outputs(
        transcripts_path=tx_out,
        boundaries_path=bd_out,
        chunksize=args.chunksize,
    )

    print("\n[done] Segger-ready files created:")
    print(f"  - {tx_out}")
    print(f"  - {bd_out}")


if __name__ == "__main__":
    main()
