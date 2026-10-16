#!/usr/bin/env python3

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


def parse_args():
    ap = argparse.ArgumentParser(
        description="Generate a 3x3 micron-to-pixel affine CSV for Xenium morphology images."
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output CSV path for the 3x3 micron-to-pixel affine matrix.",
    )
    ap.add_argument(
        "--ome-tif",
        help="Optional Xenium morphology OME-TIFF. If set, read PhysicalSizeX/Y from OME metadata.",
    )
    ap.add_argument(
        "--image-size-px",
        nargs=2,
        type=int,
        metavar=("WIDTH_PX", "HEIGHT_PX"),
        help="Optional image size in pixels. If provided, print the image bounds in microns.",
    )
    ap.add_argument(
        "--um-per-px-x",
        type=float,
        help="Microns per pixel along X. If omitted, defaults to Xenium's common value 0.2125 unless --ome-tif is used.",
    )
    ap.add_argument(
        "--um-per-px-y",
        type=float,
        help="Microns per pixel along Y. If omitted, defaults to X value.",
    )
    ap.add_argument(
        "--crop-origin-px",
        nargs=2,
        type=float,
        metavar=("X0_PX", "Y0_PX"),
        default=(0.0, 0.0),
        help="Pixel origin of the cropped ROI within the full morphology image. Use 0 0 for the full uncropped image.",
    )
    ap.add_argument(
        "--crop-origin-um",
        nargs=2,
        type=float,
        metavar=("X0_UM", "Y0_UM"),
        help="Alternative to --crop-origin-px. Useful when your crop starts at a known micron coordinate and there is no extra margin.",
    )
    ap.add_argument(
        "--pixel-bounds-file",
        help="Optional file containing text like 'px_min=..., py_min=...'. If set, this supplies the crop origin in pixels.",
    )
    return ap.parse_args()


def read_ome_pixel_size(ome_tif_path):
    try:
        from tifffile import TiffFile
    except ImportError as exc:
        raise RuntimeError("tifffile is required to read OME metadata.") from exc

    with TiffFile(ome_tif_path) as tif:
        ome_xml = tif.ome_metadata
    if not ome_xml:
        raise ValueError(f"No OME metadata found in {ome_tif_path}")

    root = ET.fromstring(ome_xml)
    pixels = root.find(".//{*}Pixels")
    if pixels is None:
        raise ValueError(f"Could not find OME Pixels metadata in {ome_tif_path}")

    sx = pixels.attrib.get("PhysicalSizeX")
    sy = pixels.attrib.get("PhysicalSizeY")
    if sx is None or sy is None:
        raise ValueError(f"PhysicalSizeX/Y missing in OME metadata for {ome_tif_path}")
    return float(sx), float(sy)


def read_ome_image_size(ome_tif_path):
    try:
        from tifffile import TiffFile
    except ImportError as exc:
        raise RuntimeError("tifffile is required to read OME metadata.") from exc

    with TiffFile(ome_tif_path) as tif:
        ome_xml = tif.ome_metadata
    if not ome_xml:
        raise ValueError(f"No OME metadata found in {ome_tif_path}")

    root = ET.fromstring(ome_xml)
    pixels = root.find(".//{*}Pixels")
    if pixels is None:
        raise ValueError(f"Could not find OME Pixels metadata in {ome_tif_path}")

    size_x = pixels.attrib.get("SizeX")
    size_y = pixels.attrib.get("SizeY")
    if size_x is None or size_y is None:
        raise ValueError(f"SizeX/SizeY missing in OME metadata for {ome_tif_path}")
    return int(size_x), int(size_y)


def read_pixel_bounds(path):
    text = Path(path).read_text().strip()
    x_match = re.search(r"px_min=([-+0-9.eE]+)", text)
    y_match = re.search(r"py_min=([-+0-9.eE]+)", text)
    if x_match is None or y_match is None:
        raise ValueError(f"Could not parse px_min/py_min from {path}")
    return float(x_match.group(1)), float(y_match.group(1))


def main():
    args = parse_args()

    if args.ome_tif:
        um_per_px_x, um_per_px_y = read_ome_pixel_size(args.ome_tif)
    else:
        um_per_px_x = args.um_per_px_x if args.um_per_px_x is not None else 0.2125
        um_per_px_y = args.um_per_px_y if args.um_per_px_y is not None else um_per_px_x

    if um_per_px_x <= 0 or um_per_px_y <= 0:
        raise ValueError("Microns per pixel must be positive.")

    crop_origin_sources = sum(
        value is not None
        for value in (args.crop_origin_um, args.pixel_bounds_file)
    )
    if crop_origin_sources > 1:
        raise ValueError(
            "Use at most one of --crop-origin-um or --pixel-bounds-file together with the default/explicit --crop-origin-px."
        )

    if args.pixel_bounds_file:
        crop_x_px, crop_y_px = read_pixel_bounds(args.pixel_bounds_file)
    elif args.crop_origin_um:
        crop_x_um, crop_y_um = args.crop_origin_um
        crop_x_px = crop_x_um / um_per_px_x
        crop_y_px = crop_y_um / um_per_px_y
    else:
        crop_x_px, crop_y_px = map(float, args.crop_origin_px)

    affine = np.array(
        [
            [1.0 / um_per_px_x, 0.0, -crop_x_px],
            [0.0, 1.0 / um_per_px_y, -crop_y_px],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out_path, affine, fmt="%.10f")

    image_size_px = None
    if args.image_size_px:
        image_size_px = tuple(args.image_size_px)
    elif args.ome_tif:
        image_size_px = read_ome_image_size(args.ome_tif)

    print(f"Wrote affine matrix to {out_path}")
    print(f"um_per_px_x={um_per_px_x:.10f}")
    print(f"um_per_px_y={um_per_px_y:.10f}")
    print(f"crop_origin_px=({crop_x_px:.10f}, {crop_y_px:.10f})")
    if image_size_px:
        width_px, height_px = image_size_px
        x_min_um = crop_x_px * um_per_px_x
        y_min_um = crop_y_px * um_per_px_y
        x_max_um = (crop_x_px + width_px) * um_per_px_x
        y_max_um = (crop_y_px + height_px) * um_per_px_y
        print(f"image_size_px=({width_px}, {height_px})")
        print(
            "roi_um_bounds="
            f"({x_min_um:.10f}, {y_min_um:.10f}, {x_max_um:.10f}, {y_max_um:.10f})"
        )
    print("Matrix:")
    np.set_printoptions(precision=10, suppress=False)
    print(affine)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
