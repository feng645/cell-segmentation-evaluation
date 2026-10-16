#!/usr/bin/env python3

import argparse
import gzip
import json
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops


def parse_args():
    ap = argparse.ArgumentParser(
        description="Prepare Cellpose outputs for the R evaluation pipeline."
    )
    ap.add_argument("--mask-gz", required=True, help="Expanded Cellpose mask as gzipped .npy")
    ap.add_argument("--transcripts", required=True, help="ROI transcript CSV")
    coord_group = ap.add_mutually_exclusive_group(required=True)
    coord_group.add_argument(
        "--roi-um",
        nargs=4,
        type=float,
        metavar=("X_MIN", "Y_MIN", "X_MAX", "Y_MAX"),
        help="ROI bounds in microns corresponding to the image."
    )
    coord_group.add_argument(
        "--affine-matrix",
        help="Optional 3x3 micron-to-pixel affine CSV for the ROI image."
    )
    ap.add_argument("--out-dir", required=True, help="Directory for Cellpose raw evaluation artifacts")
    ap.add_argument(
        "--out-label-tif",
        default=None,
        help="Optional label TIFF path. Useful for BANKSY and spatial fallbacks."
    )
    return ap.parse_args()


def load_mask(path):
    with gzip.open(path, "rb") as handle:
        mask = np.load(handle)
    if mask.ndim != 2:
        raise ValueError(f"Expected a 2D label mask, got shape {mask.shape}")
    return mask


def affine_from_roi_um(mask_shape, roi_um):
    x_min, y_min, x_max, y_max = roi_um
    height, width = mask_shape
    um_per_px_x = (x_max - x_min) / width
    um_per_px_y = (y_max - y_min) / height
    return np.array([
        [1.0 / um_per_px_x, 0.0, -x_min / um_per_px_x],
        [0.0, 1.0 / um_per_px_y, -y_min / um_per_px_y],
        [0.0, 0.0, 1.0],
    ], dtype=float)


def load_affine(path):
    affine = np.loadtxt(path)
    if affine.shape != (3, 3):
        raise ValueError(f"Expected a 3x3 affine matrix, got shape {affine.shape}")
    return affine


def resolve_affine(mask_shape, roi_um, affine_matrix):
    affine = load_affine(affine_matrix) if affine_matrix is not None else affine_from_roi_um(mask_shape, roi_um)
    inv_affine = np.linalg.inv(affine)
    return affine, inv_affine


def micron_to_pixel(affine, x_um, y_um):
    pts = np.column_stack([x_um, y_um, np.ones(len(x_um), dtype=float)])
    px = pts @ affine.T
    return px[:, 0], px[:, 1]


def pixel_to_micron(inv_affine, x_px, y_px):
    pts = np.column_stack([x_px, y_px, np.ones(len(x_px), dtype=float)])
    um = pts @ inv_affine.T
    return um[:, 0], um[:, 1]


def normalize_transcript_columns(transcripts):
    candidates = [
        ("global_x", "global_y", "gene"),
        ("x_location", "y_location", "feature_name"),
        ("x", "y", "gene"),
    ]
    for x_col, y_col, gene_col in candidates:
        if {x_col, y_col, gene_col}.issubset(transcripts.columns):
            out = transcripts.copy()
            rename_map = {}
            if x_col != "global_x":
                rename_map[x_col] = "global_x"
            if y_col != "global_y":
                rename_map[y_col] = "global_y"
            if gene_col != "gene":
                rename_map[gene_col] = "gene"
            if rename_map:
                out = out.rename(columns=rename_map)
            return out
    raise KeyError(
        "Transcript table missing a supported coordinate/gene column set. "
        "Expected one of: "
        "(global_x, global_y, gene), "
        "(x_location, y_location, feature_name), "
        "(x, y, gene)."
    )


def assign_transcripts(transcripts, mask, affine):
    transcripts = normalize_transcript_columns(transcripts)

    x_px_f, y_px_f = micron_to_pixel(
        affine,
        transcripts["global_x"].to_numpy(dtype=float),
        transcripts["global_y"].to_numpy(dtype=float),
    )
    x_px = np.floor(x_px_f).astype(np.int64)
    y_px = np.floor(y_px_f).astype(np.int64)

    height, width = mask.shape

    in_bounds = (
        (x_px >= 0) & (x_px < width) &
        (y_px >= 0) & (y_px < height)
    )
    labels = np.zeros(len(transcripts), dtype=np.int64)
    labels[in_bounds] = mask[y_px[in_bounds], x_px[in_bounds]]

    out = transcripts.copy()
    out["cell_id_cellpose"] = np.where(labels > 0, labels.astype(str), "-1")
    out["cellpose_x_px"] = np.where(in_bounds, x_px, np.nan)
    out["cellpose_y_px"] = np.where(in_bounds, y_px, np.nan)
    return out


def build_expr_matrix(assigned):
    cell_df = assigned.loc[assigned["cell_id_cellpose"] != "-1", ["cell_id_cellpose", "gene"]].copy()
    if cell_df.empty:
        return pd.DataFrame(columns=["cell_id"])
    counts = (
        cell_df
        .groupby(["cell_id_cellpose", "gene"])
        .size()
        .rename("count")
        .reset_index()
    )
    wide = counts.pivot(index="cell_id_cellpose", columns="gene", values="count").fillna(0).astype(np.int64)
    wide = wide.sort_index(axis=0).sort_index(axis=1)
    wide.index.name = "cell_id"
    return wide.reset_index()


def build_centroids(mask, inv_affine):
    rows, cols = np.nonzero(mask > 0)
    if len(rows) == 0:
        return pd.DataFrame(columns=["cell_id", "x", "y", "area_px", "area_um2"])

    work = pd.DataFrame({
        "cell_id": mask[rows, cols].astype(np.int64),
        "row": rows.astype(np.float64),
        "col": cols.astype(np.float64),
    })
    centroids = (
        work.groupby("cell_id")
        .agg(row=("row", "mean"), col=("col", "mean"), area_px=("cell_id", "size"))
        .reset_index()
    )
    x_um, y_um = pixel_to_micron(inv_affine, centroids["col"].to_numpy(), centroids["row"].to_numpy())
    um_step_x = inv_affine[:2, :2] @ np.array([1.0, 0.0])
    um_step_y = inv_affine[:2, :2] @ np.array([0.0, 1.0])
    pixel_area_um2 = float(abs(np.linalg.det(inv_affine[:2, :2])))
    centroids["x"] = x_um
    centroids["y"] = y_um
    centroids["area_um2"] = centroids["area_px"] * pixel_area_um2
    centroids["cell_id"] = centroids["cell_id"].astype(str)
    return centroids[["cell_id", "x", "y", "area_px", "area_um2"]]


def contour_to_ring(contour, inv_affine):
    pts = contour[:, 0, :]
    if len(pts) < 3:
        return None
    x_um, y_um = pixel_to_micron(inv_affine, pts[:, 0].astype(float), pts[:, 1].astype(float))
    ring = [[float(x), float(y)] for x, y in zip(x_um, y_um)]
    if ring[0] != ring[-1]:
        ring.append(ring[0])
    return ring


def build_geojson(mask, inv_affine):
    features = []
    # One regionprops pass over the image gives every label's bounding box,
    # avoiding a full-image binary array per label (previously O(n_cells x image_size))
    props = regionprops(mask)
    for prop in props:
        label = prop.label
        r0, c0, r1, c1 = prop.bbox  # bounding box in pixel coords
        # Run findContours on the bounding-box patch only, not the whole image
        patch = (mask[r0:r1, c0:c1] == label).astype(np.uint8)
        contours, _ = cv2.findContours(patch, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        rings = []
        for contour in contours:
            pts = contour[:, 0, :]
            if len(pts) < 3:
                continue
            # Add the bounding-box offset to get full-image pixel coordinates, then convert to microns
            x_px = pts[:, 0].astype(float) + c0
            y_px = pts[:, 1].astype(float) + r0
            x_um, y_um = pixel_to_micron(inv_affine, x_px, y_px)
            ring = [[float(x), float(y)] for x, y in zip(x_um, y_um)]
            if ring[0] != ring[-1]:
                ring.append(ring[0])
            rings.append(ring)

        if not rings:
            continue

        if len(rings) == 1:
            geometry = {"type": "Polygon", "coordinates": [rings[0]]}
        else:
            geometry = {"type": "MultiPolygon", "coordinates": [[ring] for ring in rings]}

        features.append({
            "type": "Feature",
            "id": str(int(label)),
            "properties": {"cell_id": str(int(label))},
            "geometry": geometry,
        })
    return {"type": "FeatureCollection", "features": features}


def main():
    args = parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    mask = load_mask(args.mask_gz)
    affine, inv_affine = resolve_affine(mask.shape, args.roi_um, args.affine_matrix)

    transcripts = pd.read_csv(args.transcripts)
    assigned = assign_transcripts(transcripts, mask, affine)
    expr_mat = build_expr_matrix(assigned)
    centroids = build_centroids(mask, inv_affine)
    geojson = build_geojson(mask, inv_affine)

    assigned_path = out_dir / "cellpose_transcript_assignments.csv"
    expr_path = out_dir / "cellpose_expr_mat.csv"
    centroids_path = out_dir / "cellpose_centroids.csv"
    geojson_path = out_dir / "cellpose_boundaries.geojson.gz"
    stats_path = out_dir / "cellpose_mask_stats.csv"

    assigned.to_csv(assigned_path, index=False)
    expr_mat.to_csv(expr_path, index=False)
    centroids.to_csv(centroids_path, index=False)
    with gzip.open(geojson_path, "wt", encoding="utf-8") as handle:
        json.dump(geojson, handle, separators=(",", ":"))

    stats = pd.DataFrame({
        "metric": ["n_cells", "n_foreground_pixels", "n_assigned_transcripts"],
        "value": [
            int(mask.max()),
            int(np.count_nonzero(mask)),
            int((assigned["cell_id_cellpose"] != "-1").sum()),
        ],
    })
    stats.to_csv(stats_path, index=False)

    if args.out_label_tif:
        tifffile.imwrite(args.out_label_tif, mask.astype(np.uint32), photometric="minisblack")

    print(f"Wrote {assigned_path}")
    print(f"Wrote {expr_path}")
    print(f"Wrote {centroids_path}")
    print(f"Wrote {geojson_path}")
    print(f"Wrote {stats_path}")


if __name__ == "__main__":
    main()
