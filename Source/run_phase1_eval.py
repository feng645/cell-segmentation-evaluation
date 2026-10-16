#!/usr/bin/env python3

import argparse
import gzip
import json
import math
import os
import re
import shutil
import sys
import gc
from collections import defaultdict

import numpy as np
import pandas as pd
from PIL import Image
from shapely import points as shapely_points
from shapely.geometry import MultiPolygon, Polygon, mapping, shape
from shapely.strtree import STRtree

Image.MAX_IMAGE_PIXELS = None


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def read_json(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return json.load(handle)
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def write_geojson_gz(path, features):
    payload = {"type": "FeatureCollection", "features": features}
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(payload, handle)


def normalize_cell_id(value):
    if value is None:
        return "-1"
    if isinstance(value, float) and math.isnan(value):
        return "-1"
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null"}:
        return "-1"
    if text in {"0", "UNASSIGNED"}:
        return "-1"
    return text


def compression_for_csv(path):
    return "gzip" if path.endswith(".gz") else None


def bool_from_value(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y"}


def feature_cell_id(feature):
    if feature.get("id") is not None:
        return str(feature["id"])
    props = feature.get("properties") or {}
    for key in ("cell_id", "cell", "id"):
        if props.get(key) is not None:
            return str(props[key])
    raise ValueError("Could not resolve cell id from feature.")


def to_polygon_geometries(geom):
    if geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return list(geom.geoms)
    return []


def load_geojson_polygons(path):
    payload = read_json(path)
    polygons = []
    features = []
    for feature in payload.get("features", []):
        cell_id = feature_cell_id(feature)
        geom = shape(feature["geometry"])
        for poly in to_polygon_geometries(geom):
            if poly.is_empty or poly.area <= 0:
                continue
            polygons.append((cell_id, poly))
        features.append({
            "type": "Feature",
            "id": cell_id,
            "properties": {"cell_id": cell_id},
            "geometry": feature["geometry"],
        })
    return polygons, features


def load_xenium_vertex_csv(path):
    frame = pd.read_csv(path)
    required = {"cell_id", "vertex_x", "vertex_y"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing Xenium boundary columns: {sorted(missing)}")

    polygons = []
    features = []
    for cell_id, sub in frame.groupby("cell_id", sort=False):
        coords = sub[["vertex_x", "vertex_y"]].to_numpy(dtype=float)
        if len(coords) < 3:
            continue
        if not np.array_equal(coords[0], coords[-1]):
            coords = np.vstack([coords, coords[0]])
        poly = Polygon(coords)
        if poly.is_empty or poly.area <= 0:
            continue
        cell_id = str(cell_id)
        polygons.append((cell_id, poly))
        features.append({
            "type": "Feature",
            "id": cell_id,
            "properties": {"cell_id": cell_id},
            "geometry": mapping(poly),
        })
    return polygons, features


def read_bidcell_affine_summary(path):
    frame = pd.read_csv(path, sep="\t", header=None, names=["key", "value"])
    values = {}
    for _, row in frame.iterrows():
        key = str(row["key"])
        raw = row["value"]
        if isinstance(raw, str) and raw.strip() in {"True", "False"}:
            values[key] = raw.strip() == "True"
        else:
            try:
                values[key] = float(raw)
            except Exception:
                values[key] = raw
    return values


def transform_bidcell_polygon(poly, transform_cfg, method_cfg):
    transform_type = ((transform_cfg or {}).get("type")) or "identity"
    if transform_type == "identity":
        return poly
    if transform_type != "bidcell_crop_tiff":
        raise ValueError(f"Unsupported BIDCell transform type: {transform_type}")

    affine_summary = method_cfg["input"].get("affine_summary")
    seg_path = method_cfg["input"].get("seg_path")
    if not affine_summary or not seg_path:
        raise ValueError("BIDCell crop transform requires affine_summary and seg_path.")

    affine = read_bidcell_affine_summary(affine_summary)
    with Image.open(seg_path) as seg_img:
        width, height = seg_img.size

    x_scale = float(affine["size_x"]) / float(width)
    y_scale = float(affine["size_y"]) / float(height)
    min_x = float(affine["min_x"]) + float(affine.get("global_shift_x", 0.0))
    min_y = float(affine["min_y"]) + float(affine.get("global_shift_y", 0.0))
    ts_x = float(affine.get("scale_ts_x", 1.0))
    ts_y = float(affine.get("scale_ts_y", 1.0))

    def transform_ring(coords):
        arr = np.asarray(coords, dtype=float)
        arr[:, 0] = min_x + arr[:, 0] * x_scale * ts_x
        arr[:, 1] = min_y + arr[:, 1] * y_scale * ts_y
        return arr

    exterior = transform_ring(poly.exterior.coords)
    interiors = [transform_ring(ring.coords) for ring in poly.interiors]
    return Polygon(exterior, interiors)


def load_bidcell_geojson(path, method_cfg):
    payload = read_json(path)
    polygons = []
    features = []
    transform_cfg = method_cfg["input"].get("transform") or {}
    for feature in payload.get("features", []):
        cell_id = feature_cell_id(feature)
        geom = shape(feature["geometry"])
        transformed_parts = []
        for poly in to_polygon_geometries(geom):
            transformed = transform_bidcell_polygon(poly, transform_cfg, method_cfg)
            if transformed.is_empty or transformed.area <= 0:
                continue
            polygons.append((cell_id, transformed))
            transformed_parts.append(transformed)
        if not transformed_parts:
            continue
        merged_geom = transformed_parts[0] if len(transformed_parts) == 1 else MultiPolygon(transformed_parts)
        features.append({
            "type": "Feature",
            "id": cell_id,
            "properties": {"cell_id": cell_id},
            "geometry": mapping(merged_geom),
        })
    return polygons, features


def prepare_polygon_method(method_cfg, polygons_dir):
    loader = method_cfg["input"]["loader"]
    path = method_cfg["input"].get("path")
    if not path or not os.path.exists(path):
        return {"name": method_cfg["name"], "status": "missing_input", "reason": f"Missing input path: {path}"}

    if loader == "geojson":
        polygons, features = load_geojson_polygons(path)
    elif loader == "xenium_vertex_csv":
        polygons, features = load_xenium_vertex_csv(path)
    elif loader == "bidcell_geojson":
        polygons, features = load_bidcell_geojson(path, method_cfg)
    else:
        raise ValueError(f"Unsupported polygon loader: {loader}")

    if not polygons:
        return {"name": method_cfg["name"], "status": "empty", "reason": "No usable polygons were loaded."}

    polygon_out = os.path.join(polygons_dir, f"{method_cfg['name']}.geojson.gz")
    write_geojson_gz(polygon_out, features)

    geometries = [poly for _, poly in polygons]
    tree = STRtree(geometries)
    ids = np.asarray([cell_id for cell_id, _ in polygons], dtype=object)
    areas = np.asarray([poly.area for _, poly in polygons], dtype=float)
    area_by_cell = defaultdict(float)
    for cell_id, poly in polygons:
        area_by_cell[cell_id] += float(poly.area)

    return {
        "name": method_cfg["name"],
        "mode": "polygon",
        "status": "ready",
        "loader": loader,
        "polygon_path": polygon_out,
        "tree": tree,
        "cell_ids": ids,
        "areas": areas,
        "cell_areas": dict(area_by_cell),
        "n_cells_raw": int(len(area_by_cell)),
    }


def load_mask_method(method_cfg):
    path = method_cfg["input"].get("path")
    affine_path = method_cfg["input"].get("affine_matrix")
    if not path or not os.path.exists(path):
        return {"name": method_cfg["name"], "status": "missing_input", "reason": f"Missing mask path: {path}"}
    if not affine_path or not os.path.exists(affine_path):
        return {"name": method_cfg["name"], "status": "missing_input", "reason": f"Missing affine matrix path: {affine_path}"}

    with Image.open(path) as img:
        mask = np.array(img)
    affine = pd.read_csv(affine_path, header=None, sep=r"[\s,]+", engine="python").to_numpy(dtype=float)
    if affine.shape != (3, 3):
        raise ValueError(f"Affine matrix must be 3x3: {affine_path}")
    unique_labels = np.unique(mask)
    unique_labels = unique_labels[np.rint(unique_labels).astype(np.int64) > 0]

    return {
        "name": method_cfg["name"],
        "mode": "mask",
        "status": "ready",
        "loader": method_cfg["input"]["loader"],
        "mask": mask,
        "affine": affine,
        "height": int(mask.shape[0]),
        "width": int(mask.shape[1]),
        "n_cells_raw": int(len(unique_labels)),
    }


def prepare_method(method_cfg, polygons_dir):
    if not method_cfg.get("enabled", True):
        return {
            "name": method_cfg["name"],
            "status": "disabled",
            "reason": method_cfg.get("reason", "Disabled by run config."),
        }

    loader = method_cfg["input"]["loader"]
    if loader in {"geojson", "xenium_vertex_csv", "bidcell_geojson"}:
        return prepare_polygon_method(method_cfg, polygons_dir)
    if loader == "mask_tiff_affine":
        return load_mask_method(method_cfg)
    raise ValueError(f"Unsupported method loader: {loader}")


def standardize_transcript_chunk(chunk, column_cfg, state):
    out = pd.DataFrame()

    transcript_col = column_cfg.get("transcript_id")
    if transcript_col and transcript_col in chunk.columns:
        out["transcript_id"] = chunk[transcript_col].astype(str)
    else:
        start = state["generated_transcript_id"]
        stop = start + len(chunk)
        out["transcript_id"] = [f"generated_{idx}" for idx in range(start, stop)]
        state["generated_transcript_id"] = stop

    def optional_column(col_name):
        if col_name and col_name in chunk.columns:
            return chunk[col_name]
        return pd.Series([np.nan] * len(chunk), index=chunk.index)

    def optional_text_column(col_name):
        series = optional_column(col_name)
        return series.where(series.notna(), "")

    required_cols = {
        "gene": column_cfg.get("gene"),
        "x": column_cfg.get("x"),
        "y": column_cfg.get("y"),
    }
    for canonical, src in required_cols.items():
        if not src or src not in chunk.columns:
            raise ValueError(f"Missing required transcript column mapping for `{canonical}` -> `{src}`")

    out["gene"] = chunk[column_cfg["gene"]].astype(str)
    out["x"] = pd.to_numeric(chunk[column_cfg["x"]], errors="coerce")
    out["y"] = pd.to_numeric(chunk[column_cfg["y"]], errors="coerce")
    out["z"] = pd.to_numeric(optional_column(column_cfg.get("z")), errors="coerce")
    out["prior_cell_id"] = optional_column(column_cfg.get("prior_cell_id")).map(normalize_cell_id)
    out["fov"] = optional_text_column(column_cfg.get("fov")).astype(str)
    out["qv"] = pd.to_numeric(optional_column(column_cfg.get("qv")), errors="coerce")
    out["overlaps_nucleus"] = optional_column(column_cfg.get("overlaps_nucleus")).map(bool_from_value)

    return out


def apply_transcript_filters(chunk, filter_cfg):
    keep = np.ones(len(chunk), dtype=bool)

    if filter_cfg.get("drop_blank_genes", True):
        keep &= chunk["gene"].notna().to_numpy()
        keep &= (chunk["gene"].astype(str).str.strip() != "").to_numpy()

    min_qv = filter_cfg.get("min_qv")
    if min_qv is not None and "qv" in chunk.columns:
        keep &= np.nan_to_num(chunk["qv"].to_numpy(), nan=-np.inf) >= float(min_qv)

    exclude_regex = filter_cfg.get("exclude_gene_regex")
    if exclude_regex:
        pattern = re.compile(exclude_regex)
        keep &= ~chunk["gene"].astype(str).map(lambda x: bool(pattern.search(x))).to_numpy()

    keep &= np.isfinite(chunk["x"].to_numpy()) & np.isfinite(chunk["y"].to_numpy())
    return chunk.loc[keep].copy()


def standardize_transcripts(config, transcripts_dir):
    input_path = config["canonical_transcripts"]["path"]
    chunk_size = int(config["runtime"]["chunk_size"])
    column_cfg = config["platform_protocol"]["transcript_columns"]
    filter_cfg = config["platform_protocol"].get("transcript_filters") or {}

    source_path = os.path.join(transcripts_dir, "transcripts_source.csv.gz")
    eval_path = os.path.join(transcripts_dir, "transcripts_eval.csv.gz")

    source_header_written = False
    eval_header_written = False

    state = {"generated_transcript_id": 1}
    n_source = 0
    n_eval = 0

    read_kwargs = {"chunksize": chunk_size}
    if input_path.endswith(".gz"):
        read_kwargs["compression"] = "gzip"

    print(f"[phase1] standardizing transcripts from {input_path}", flush=True)
    with gzip.open(source_path, "wt", encoding="utf-8", newline="") as source_writer, \
         gzip.open(eval_path, "wt", encoding="utf-8", newline="") as eval_writer:
        for chunk_idx, chunk in enumerate(pd.read_csv(input_path, **read_kwargs), start=1):
            standardized = standardize_transcript_chunk(chunk, column_cfg, state)
            filtered = apply_transcript_filters(standardized, filter_cfg)

            standardized.to_csv(source_writer, index=False, header=not source_header_written)
            filtered.to_csv(eval_writer, index=False, header=not eval_header_written)
            source_header_written = True
            eval_header_written = True
            n_source += int(len(standardized))
            n_eval += int(len(filtered))

            if chunk_idx == 1 or chunk_idx % 10 == 0:
                print(
                    f"[phase1] transcript chunks={chunk_idx} source_rows={n_source} eval_rows={n_eval}",
                    flush=True,
                )

    return {
        "source_path": source_path,
        "eval_path": eval_path,
        "n_source_transcripts": n_source,
        "n_eval_transcripts": n_eval,
    }


def assign_chunk_polygon(method_state, frame):
    assignments = np.full(len(frame), "-1", dtype=object)
    if len(frame) == 0:
        return assignments

    pts = shapely_points(frame["x"].to_numpy(dtype=float), frame["y"].to_numpy(dtype=float))
    left_idx, right_idx = method_state["tree"].query(pts, predicate="intersects")
    if len(left_idx) == 0:
        return assignments

    poly_areas = method_state["areas"][right_idx]
    order = np.lexsort((poly_areas, left_idx))
    left_sorted = left_idx[order]
    right_sorted = right_idx[order]
    first = np.ones(len(left_sorted), dtype=bool)
    first[1:] = left_sorted[1:] != left_sorted[:-1]
    assignments[left_sorted[first]] = method_state["cell_ids"][right_sorted[first]]
    return assignments


def assign_chunk_mask(method_state, frame):
    assignments = np.full(len(frame), "-1", dtype=object)
    if len(frame) == 0:
        return assignments

    pts = np.column_stack([
        frame["x"].to_numpy(dtype=float),
        frame["y"].to_numpy(dtype=float),
        np.ones(len(frame), dtype=float),
    ])
    pixel = pts @ method_state["affine"].T
    px = np.floor(pixel[:, 0]).astype(int)
    py = np.floor(pixel[:, 1]).astype(int)

    in_bounds = (px >= 0) & (px < method_state["width"]) & (py >= 0) & (py < method_state["height"])
    if not np.any(in_bounds):
        return assignments

    labels = method_state["mask"][py[in_bounds], px[in_bounds]]
    labels = np.rint(labels).astype(np.int64)
    label_text = labels.astype(str)
    label_text[labels <= 0] = "-1"
    assignments[in_bounds] = label_text
    return assignments


def lightweight_method_record(method_state):
    return {
        "name": method_state["name"],
        "status": method_state.get("status"),
        "reason": method_state.get("reason"),
        "mode": method_state.get("mode"),
        "loader": method_state.get("loader"),
        "polygon_path": method_state.get("polygon_path"),
        "n_cells_raw": method_state.get("n_cells_raw"),
        "cell_areas": method_state.get("cell_areas", {}),
    }


def append_method_assignments(assignments_path, method_state, chunk_size):
    column_name = f"cell_id_{method_state['name']}"
    temp_path = assignments_path + ".tmp"
    summary = {
        "method": method_state["name"],
        "status": method_state.get("status"),
        "n_cells_raw": int(method_state.get("n_cells_raw", 0) or 0),
        "assigned_transcripts": 0,
        "background_transcripts": 0,
        "polygon_path": method_state.get("polygon_path"),
        "mode": method_state.get("mode"),
        "reason": method_state.get("reason"),
    }

    header_written = False
    try:
        with gzip.open(temp_path, "wt", encoding="utf-8", newline="") as writer:
            for chunk_idx, chunk in enumerate(pd.read_csv(assignments_path, chunksize=chunk_size, compression="gzip"), start=1):
                chunk_out = chunk.copy()
                if method_state.get("status") == "ready":
                    if method_state["mode"] == "polygon":
                        assignments = assign_chunk_polygon(method_state, chunk_out)
                    elif method_state["mode"] == "mask":
                        assignments = assign_chunk_mask(method_state, chunk_out)
                    else:
                        raise ValueError(f"Unsupported method mode: {method_state['mode']}")
                    assigned = assignments != "-1"
                    summary["assigned_transcripts"] += int(np.sum(assigned))
                    summary["background_transcripts"] += int(np.sum(~assigned))
                else:
                    assignments = np.full(len(chunk_out), "-1", dtype=object)
                    summary["background_transcripts"] += int(len(chunk_out))

                chunk_out[column_name] = assignments
                chunk_out.to_csv(writer, index=False, header=not header_written)
                header_written = True

                if chunk_idx == 1 or chunk_idx % 10 == 0:
                    print(
                        f"[phase1] method={method_state['name']} chunks={chunk_idx} assigned={summary['assigned_transcripts']}",
                        flush=True,
                    )

        os.replace(temp_path, assignments_path)
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)
    return summary


def build_assignments(config, transcript_meta, assignments_dir, polygons_dir):
    eval_path = transcript_meta["eval_path"]
    output_path = os.path.join(assignments_dir, "transcripts_with_assignments.csv.gz")
    summary_path = os.path.join(assignments_dir, "assignment_summary_by_method.csv")
    chunk_size = int(config["runtime"]["chunk_size"])

    shutil.copyfile(eval_path, output_path)
    summaries = []
    method_records = []

    for method_cfg in config["methods"]:
        print(f"[phase1] preparing method={method_cfg['name']}", flush=True)
        method_state = prepare_method(method_cfg, polygons_dir)
        summary = append_method_assignments(output_path, method_state, chunk_size)
        print(
            f"[phase1] completed method={summary['method']} status={summary['status']} assigned={summary['assigned_transcripts']}",
            flush=True,
        )
        summaries.append(summary)
        method_records.append(lightweight_method_record(method_state))
        del method_state
        gc.collect()

    pd.DataFrame(summaries).sort_values("method").to_csv(summary_path, index=False)
    return {
        "assignments_path": output_path,
        "summary_path": summary_path,
        "summary": summaries,
        "method_records": method_records,
    }


def aggregate_method_counts(assignments_path, method_name, chunk_size):
    assignment_col = f"cell_id_{method_name}"
    counts = defaultdict(int)

    for chunk in pd.read_csv(assignments_path, chunksize=chunk_size, compression="gzip"):
        if assignment_col not in chunk.columns:
            raise ValueError(f"Assignment column not found: {assignment_col}")
        sub = chunk.loc[chunk[assignment_col].astype(str) != "-1", [assignment_col, "gene"]].copy()
        if sub.empty:
            continue
        grouped = sub.groupby([assignment_col, "gene"]).size()
        for (cell_id, gene), count in grouped.items():
            counts[(str(cell_id), str(gene))] += int(count)

    return counts


def build_cell_outputs(config, assignments_meta, methods_state, expression_dir, cell_tables_dir):
    chunk_size = int(config["runtime"]["chunk_size"])
    mt_regex = re.compile(config["platform_protocol"].get("mt_gene_regex") or "(?i)^mt[-_]")
    cell_filters = config["platform_protocol"].get("cell_filters") or {}
    min_genes = int(cell_filters.get("min_genes", 30))
    max_genes = int(cell_filters.get("max_genes", 8000))
    max_pct_mt = float(cell_filters.get("max_pct_mt", 15.0))

    manifest_rows = []

    for method in methods_state:
        method_name = method["name"]
        if method.get("status") != "ready":
            manifest_rows.append({
                "method": method_name,
                "status": method.get("status"),
                "reason": method.get("reason"),
                "n_cells_raw": 0,
                "n_cells_eval": 0,
                "counts_long_path": None,
                "cell_table_path": None,
            })
            continue

        print(f"[phase1] aggregating cell outputs for method={method_name}", flush=True)
        counts = aggregate_method_counts(assignments_meta["assignments_path"], method_name, chunk_size)
        counts_rows = [
            {"cell_id": cell_id, "gene": gene, "count": count}
            for (cell_id, gene), count in counts.items()
        ]
        counts_df = pd.DataFrame(counts_rows)

        counts_long_path = os.path.join(expression_dir, f"{method_name}_counts_long.csv.gz")
        if counts_df.empty:
            counts_df = pd.DataFrame(columns=["cell_id", "gene", "count"])
        counts_df.to_csv(counts_long_path, index=False, compression="gzip")

        if counts_df.empty:
            cell_df = pd.DataFrame(columns=["cell_id", "n_transcripts", "n_genes", "n_mt_transcripts", "pct_mt", "cell_area", "keep_for_eval"])
        else:
            counts_df["is_mt"] = counts_df["gene"].astype(str).map(lambda gene: bool(mt_regex.search(gene)))
            per_cell = counts_df.groupby("cell_id").agg(
                n_transcripts=("count", "sum"),
                n_genes=("gene", "size"),
            ).reset_index()
            mt_per_cell = counts_df.loc[counts_df["is_mt"]].groupby("cell_id")["count"].sum().rename("n_mt_transcripts")
            cell_df = per_cell.merge(mt_per_cell, on="cell_id", how="left")
            cell_df["n_mt_transcripts"] = cell_df["n_mt_transcripts"].fillna(0).astype(int)
            cell_df["pct_mt"] = np.where(
                cell_df["n_transcripts"] > 0,
                100.0 * cell_df["n_mt_transcripts"] / cell_df["n_transcripts"],
                np.nan,
            )
            if method["mode"] == "polygon":
                cell_df["cell_area"] = cell_df["cell_id"].map(method.get("cell_areas", {}))
            else:
                cell_df["cell_area"] = np.nan
            cell_df["keep_for_eval"] = (
                (cell_df["n_genes"] >= min_genes)
                & (cell_df["n_genes"] <= max_genes)
                & (cell_df["pct_mt"].fillna(0) <= max_pct_mt)
            )

        cell_table_path = os.path.join(cell_tables_dir, f"{method_name}_cells_raw.csv.gz")
        cell_df.to_csv(cell_table_path, index=False, compression="gzip")

        eval_table_path = os.path.join(cell_tables_dir, f"{method_name}_cells_eval.csv.gz")
        cell_df.loc[cell_df["keep_for_eval"] == True].to_csv(eval_table_path, index=False, compression="gzip")

        manifest_rows.append({
            "method": method_name,
            "status": "ready",
            "reason": None,
            "n_cells_raw": int(len(cell_df)),
            "n_cells_eval": int(cell_df["keep_for_eval"].fillna(False).sum()),
            "counts_long_path": counts_long_path,
            "cell_table_path": cell_table_path,
        })

    manifest_path = os.path.join(cell_tables_dir, "method_cell_output_manifest.csv")
    pd.DataFrame(manifest_rows).sort_values("method").to_csv(manifest_path, index=False)
    return manifest_path


def strip_runtime_objects(methods_state):
    stripped = []
    for method in methods_state:
        entry = {
            "name": method["name"],
            "status": method.get("status"),
            "reason": method.get("reason"),
            "mode": method.get("mode"),
            "loader": method.get("loader"),
            "polygon_path": method.get("polygon_path"),
            "n_cells_raw": method.get("n_cells_raw"),
        }
        stripped.append(entry)
    return stripped


def run_phase1(config):
    output_root = config["output_root"]
    ensure_dir(output_root)

    transcripts_dir = os.path.join(output_root, "transcripts")
    polygons_dir = os.path.join(output_root, "polygons")
    assignments_dir = os.path.join(output_root, "assignments")
    expression_dir = os.path.join(output_root, "expression")
    cell_tables_dir = os.path.join(output_root, "cell_tables")
    for path in (transcripts_dir, polygons_dir, assignments_dir, expression_dir, cell_tables_dir):
        ensure_dir(path)

    transcript_meta = standardize_transcripts(config, transcripts_dir)
    assignments_meta = build_assignments(config, transcript_meta, assignments_dir, polygons_dir)
    methods_state = assignments_meta["method_records"]
    cell_manifest_path = build_cell_outputs(config, assignments_meta, methods_state, expression_dir, cell_tables_dir)

    manifest = {
        "dataset_id": config["dataset_id"],
        "platform": config["platform"],
        "run_config_path": config.get("run_config_path"),
        "profile_path": config.get("profile_path"),
        "canonical_transcripts": {
            "input_path": config["canonical_transcripts"]["path"],
            "transcripts_source_path": transcript_meta["source_path"],
            "transcripts_eval_path": transcript_meta["eval_path"],
            "n_source_transcripts": transcript_meta["n_source_transcripts"],
            "n_eval_transcripts": transcript_meta["n_eval_transcripts"],
        },
        "assignments": {
            "transcripts_with_assignments": assignments_meta["assignments_path"],
            "assignment_summary": assignments_meta["summary_path"],
        },
        "cell_output_manifest": cell_manifest_path,
        "methods": strip_runtime_objects(methods_state),
    }
    write_json(os.path.join(output_root, "manifest.json"), manifest)


def main():
    parser = argparse.ArgumentParser(description="Run MainV2 phase1 evaluation prep.")
    parser.add_argument("--config-json", required=True, help="Compiled JSON config path.")
    args = parser.parse_args()

    with open(args.config_json, "r", encoding="utf-8") as handle:
        config = json.load(handle)

    run_phase1(config)


if __name__ == "__main__":
    main()
