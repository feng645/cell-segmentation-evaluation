#!/usr/bin/env python3

import argparse
import gzip

import numpy as np
from scipy.ndimage import distance_transform_edt
import tifffile


def parse_args():
    ap = argparse.ArgumentParser(
        description="Expand Cellpose nuclei labels by a fixed physical distance using ROI micron bounds."
    )
    ap.add_argument("--mask-gz", required=True, help="Input Cellpose mask as gzipped .npy")
    coord_group = ap.add_mutually_exclusive_group(required=True)
    coord_group.add_argument(
        "--roi-um",
        nargs=4,
        type=float,
        metavar=("X_MIN", "Y_MIN", "X_MAX", "Y_MAX"),
        help="ROI bounds in microns for the mask image."
    )
    coord_group.add_argument(
        "--affine-matrix",
        help="Optional 3x3 micron-to-pixel affine CSV for the ROI image."
    )
    ap.add_argument("--expand-um", type=float, default=5.0, help="Expansion distance in microns.")
    ap.add_argument("--out-mask-gz", required=True, help="Output expanded mask as gzipped .npy")
    ap.add_argument("--out-label-tif", default=None, help="Optional output label TIFF for downstream tools.")
    return ap.parse_args()


def load_mask(path):
    with gzip.open(path, "rb") as handle:
        mask = np.load(handle)
    if mask.ndim != 2:
        raise ValueError(f"Expected a 2D label mask, got shape {mask.shape}")
    return mask


def save_mask(path, mask):
    with gzip.open(path, "wb") as handle:
        np.save(handle, mask)


def load_affine(path):
    affine = np.loadtxt(path)
    if affine.shape != (3, 3):
        raise ValueError(f"Expected a 3x3 affine matrix, got shape {affine.shape}")
    return affine


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


def resolve_affine(mask_shape, roi_um, affine_matrix):
    if affine_matrix is not None:
        return load_affine(affine_matrix)
    return affine_from_roi_um(mask_shape, roi_um)


def main():
    args = parse_args()
    mask = load_mask(args.mask_gz)
    affine = resolve_affine(mask.shape, args.roi_um, args.affine_matrix)
    inv_affine = np.linalg.inv(affine)

    um_step_x = inv_affine[:2, :2] @ np.array([1.0, 0.0])
    um_step_y = inv_affine[:2, :2] @ np.array([0.0, 1.0])
    um_per_px_x = float(np.linalg.norm(um_step_x))
    um_per_px_y = float(np.linalg.norm(um_step_y))

    foreground = mask > 0
    if not np.any(foreground):
        raise ValueError("Input mask has no labeled pixels.")

    dist_um, nearest_idx = distance_transform_edt(
        ~foreground,
        sampling=(um_per_px_y, um_per_px_x),
        return_indices=True
    )

    expanded = mask.copy()
    grow_idx = (~foreground) & (dist_um <= args.expand_um)
    nearest_y = nearest_idx[0][grow_idx]
    nearest_x = nearest_idx[1][grow_idx]
    expanded[grow_idx] = mask[nearest_y, nearest_x]

    save_mask(args.out_mask_gz, expanded.astype(mask.dtype, copy=False))
    if args.out_label_tif:
        tifffile.imwrite(args.out_label_tif, expanded.astype(np.uint32), photometric="minisblack")

    orig_cells = int(mask.max())
    expanded_area = int(np.count_nonzero(expanded))
    orig_area = int(np.count_nonzero(mask))
    print(
        f"Expanded {orig_cells} labels by {args.expand_um:.2f} um "
        f"(pixel size: {um_per_px_x:.4f} x {um_per_px_y:.4f} um/px)."
    )
    print(f"Foreground pixels: {orig_area} -> {expanded_area}")


if __name__ == "__main__":
    main()
