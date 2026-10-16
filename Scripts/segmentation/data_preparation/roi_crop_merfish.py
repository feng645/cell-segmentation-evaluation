#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, os, re
import numpy as np
from tifffile import imread, imsave

# only needed for visualization or transcripts
import matplotlib.pyplot as plt
import pandas as pd

def parse_args():
    ap = argparse.ArgumentParser(
        description="Crop a single DAPI mosaic by a micron-space ROI; "
                    "optionally filter ROI transcripts and visualize."
    )
    ap.add_argument("--dapi", required=True, help="Path to ONE DAPI mosaic TIFF (e.g. mosaic_DAPI_z3.tif)")
    ap.add_argument("--transform", required=True, help="Path to images/micron_to_mosaic_pixel_transform.csv")
    ap.add_argument("--roi-um", nargs=4, type=float, metavar=("X_MIN","Y_MIN","X_MAX","Y_MAX"),
                    required=True, help="Rectangle ROI in MICRONS (global coordinates).")
    ap.add_argument("--out-img", required=True, help="Output cropped DAPI TIFF path.")
    ap.add_argument("--margin-px", type=int, default=64, help="Extra pixel margin (default: 64).")

    ap.add_argument("--transcripts", default=None, help="Optional detected_transcripts.csv")
    ap.add_argument("--z", type=int, default=None, help="Optional z index (override file name).")
    ap.add_argument("--viz", choices=["True","False"], default=False,
                    help="If set, export a PNG visualization for DAPI and transcripts.")
    ap.add_argument("--out-dir", default=None, help="Optional folder for aux outputs (PNG/CSV).")
    return ap.parse_args()

def load_affine_csv(path):
    A = np.loadtxt(path)
    if A.shape != (3,3):
        raise ValueError("Transform must be a 3x3 matrix.")
    return A

def um_to_px(A, x_um, y_um):
    px, py, _ = A @ np.array([x_um, y_um, 1.0], dtype=float)
    return float(px), float(py)

def infer_z_from_name(dapi_path):
    m = re.search(r"_z(\d+)\.tif$", os.path.basename(dapi_path))
    return int(m.group(1)) if m else None

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.out_img) or ".", exist_ok=True)
    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)

    # Read the image and the 3x3 affine
    img = imread(args.dapi)                   # 2D, uint16
    H, W = img.shape
    A = load_affine_csv(args.transform)       # microns -> pixels

    # Pixel crop window (affine-transform the ROI corners)
    x0, y0, x1, y1 = args.roi_um
    corners = [(x0,y0),(x0,y1),(x1,y0),(x1,y1)]
    pxs, pys = zip(*[um_to_px(A,x,y) for x,y in corners])
    xmin = max(0, int(np.floor(min(pxs))) - args.margin_px)
    xmax = min(W, int(np.ceil (max(pxs))) + args.margin_px)
    ymin = max(0, int(np.floor(min(pys))) - args.margin_px)
    ymax = min(H, int(np.ceil (max(pys))) + args.margin_px)
    if xmin >= xmax or ymin >= ymax:
        raise ValueError("ROI is outside mosaic after transform.")

    # Crop and save DAPI
    crop = img[ymin:ymax, xmin:xmax]
    imsave(args.out_img, crop)

    # Write the affine of the cropped image (for later alignment)
    A2 = A.copy()
    A2[0,2] = A[0,2] - xmin
    A2[1,2] = A[1,2] - ymin
    xform_out = (
        os.path.join(os.path.dirname(args.out_img) or ".", "roi_micron_to_pixel_transform.csv")
        if args.out_dir is None else
        os.path.join(args.out_dir, "roi_micron_to_pixel_transform.csv")
    )
    np.savetxt(xform_out, A2, fmt="%.10f")

    # Base output directory for visualizations
    viz_dir = args.out_dir or os.path.dirname(args.out_img) or "."

    # Visualize DAPI only
    if args.viz == "True":
        plt.figure(figsize=(6,6), dpi=200)
        plt.imshow(crop, cmap="gray")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(viz_dir, "viz_dapi.png"), bbox_inches="tight", pad_inches=0)
        plt.close()

    # Optional transcripts
    if args.transcripts:
        # Keep transcripts inside the ROI (micron coordinates)
        usecols = ["barcode_id","global_x","global_y","global_z","x","y","fov","gene","transcript_id","cell_id"]
        df = pd.read_csv(args.transcripts, usecols=usecols)
        m_roi = (df.global_x.between(x0, x1)) & (df.global_y.between(y0, y1))
        sub = df.loc[m_roi].copy()
        if args.out_dir:
            sub.to_csv(os.path.join(args.out_dir, "roi_transcripts.csv"), index=False)

        # Plot transcripts (pixel scatter aligned to the crop)
        if args.viz == "True" and len(sub):
            plt.figure(figsize=(6,6), dpi=200)
            plt.scatter(sub.global_x, sub.global_y, s=1)
            plt.xlim([x0, x1])
            plt.ylim([y0, y1])  # mind the y-axis direction
            plt.xlabel('Global X (μm)')
            plt.ylabel('Global Y (μm)')
            plt.title('Transcripts in ROI')
            plt.tight_layout()
            plt.savefig(os.path.join(viz_dir, "viz_transcripts.png"),
                        bbox_inches="tight", pad_inches=0)
            plt.close()

    # Record the pixel window
    if args.out_dir:
        with open(os.path.join(args.out_dir, "roi_pixel_bounds.txt"), "w") as f:
            f.write(f"px_min={xmin}, py_min={ymin}, px_max={xmax}, py_max={ymax}\n")

    print("Done.")

if __name__ == "__main__":
    main()



# Example Usage
# # Crop DAPI only (no transcripts)
# python roi_crop_merfish.py \
#   --dapi /scratch.global/feng0421/data/MERFISH/dataset/HumanLungCancerPatient1_images_mosaic_DAPI_z0.tif \
#   --transform /scratch.global/feng0421/data/MERFISH/dataset/HumanLungCancerPatient1_images_micron_to_mosaic_pixel_transform.csv \
#   --roi-um 800 800 900 900 \
#   --out-img /scratch.global/feng0421/data/MERFISH/dataset/merfish_roi/roi_dapi.tif \
#   --out-dir /scratch.global/feng0421/data/MERFISH/dataset/merfish_roi/ \
#   --viz True

# # Crop DAPI + filter ROI transcripts + DAPI PNG
# python crop_one_dapi.py \
#   --dapi /path/region_0/images/mosaic_DAPI_z3.tif \
#   --transform /path/region_0/images/micron_to_mosaic_pixel_transform.csv \
#   --transcripts /path/region_0/detected_transcripts.csv \
#   --roi-um 2000 1500 3500 2800 \
#   --out-img ./roi_dapi.tif \
#   --viz dapi

# # Visualize transcripts only (z taken from the file name, e.g. z=3; or set --z 2)
# python crop_one_dapi.py \
#   --dapi /path/region_0/images/mosaic_DAPI_z3.tif \
#   --transform /path/region_0/images/micron_to_mosaic_pixel_transform.csv \
#   --transcripts /path/region_0/detected_transcripts.csv \
#   --roi-um 2000 1500 3500 2800 \
#   --viz transcripts \
#   --out-img ./roi_dapi.tif
