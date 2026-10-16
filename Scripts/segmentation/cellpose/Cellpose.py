#!/usr/bin/env python
"""
Improved Cellpose script specifically designed for DAPI-stained images 
from Xenium, CosMx, and Merfish spatial transcriptomics platforms
"""
from argparse import ArgumentParser
from cellpose import models, io
from cellpose.io import imread
import numpy as np
import gzip
from cellpose import plot
import matplotlib.pyplot as plt
import torch
import gc

def divceil(n: int, d: int):
    """Ceiling division: returns the smallest integer >= n/d"""
    q, r = divmod(n, d)
    return q + bool(r)

def get_platform_params(platform):
    """Return optimized parameters for different spatial transcriptomics platforms"""
    params = {
        'xenium': {
            'model_type': 'nuclei',
            'diameter': None,  # Let Cellpose auto-estimate
            'channels': [0, 0],  # DAPI is single-channel
            'flow_threshold': 0.6,
            'cellprob_threshold': -2.0
        },
        'cosmx': {
            'model_type': 'nuclei', 
            'diameter': None,
            'channels': [0, 0],
            'flow_threshold': 0.8,
            'cellprob_threshold': -1.5
        },
        'merfish': {
            'model_type': 'nuclei',
            'diameter': None,
            'channels': [0, 0], 
            'flow_threshold': 0.7,
            'cellprob_threshold': -2.0
        }
    }
    return params.get(platform.lower(), params['xenium'])

if __name__ == "__main__":
    ap = ArgumentParser("Cellpose segmentation tool for DAPI-stained images from spatial transcriptomics platforms")
    ap.add_argument("input_filename", help="TIFF format DAPI-stained image")
    ap.add_argument("--platform", choices=['xenium', 'cosmx', 'merfish'], 
                   default='xenium', help="Spatial transcriptomics platform type")
    ap.add_argument("--masks-output", default="cellpose-masks.npy.gz")
    ap.add_argument("--cellprobs-output", default="cellpose-cellprobs.npy.gz")
    ap.add_argument("--tile-size", type=int, default=None, 
                   help="Split the image into squares of this size to save memory (introduces edge effects)")
    ap.add_argument("--diameter", type=int, default=None,
                   help="Manually specify nucleus diameter (overrides platform default)")
    ap.add_argument("--plot", action='store_true', help="Plot intermediate results")

    args = ap.parse_args()
    io.logger_setup()
    
    # Get platform-specific parameters
    params = get_platform_params(args.platform)
    if args.diameter is not None:
        params['diameter'] = args.diameter
    
    print(f"Using {args.platform} platform parameters: diameter={'auto-estimate' if params['diameter'] is None else params['diameter']}, "
          f"flow_threshold={params['flow_threshold']}")
    
    # Initialize model - use nuclei model for DAPI instead of default cytoplasm model
    model = models.CellposeModel(model_type=params['model_type'], gpu=True)
    
    # Read image
    img = imread(args.input_filename)
    if img.ndim == 2:
        img = np.expand_dims(img, 0)
    assert img.ndim == 3, f"Expected img to have 3 dimensions, but got {img.ndim} dimensions"
    c, w, h = img.shape
    
    # Set tile size - if None, use full image size
    if args.tile_size is None:
        args.tile_size = max(w, h)
    
    # Initialize output arrays
    masks = np.zeros((w, h), dtype=np.uint32)
    cellprob = np.zeros((w, h), dtype=np.float32)

    # Calculate number of tiles needed
    x_tiles = divceil(w, args.tile_size)
    y_tiles = divceil(h, args.tile_size)
    
    nextcell = 0  # Track next available cell ID
    
    for i in range(x_tiles):
        for j in range(y_tiles):
            print(f"Tile {i+1}/{x_tiles}, {j+1}/{y_tiles}")
            
            # Calculate tile boundaries
            x_from = i * args.tile_size
            x_to = min((i+1) * args.tile_size, w)
            y_from = j * args.tile_size
            y_to = min((j+1) * args.tile_size, h)
            
            # Extract tile from image
            tile = img[:, x_from:x_to, y_from:y_to]
            
            # Run Cellpose segmentation on the tile
            masks_ij, flows_ij, styles_ij = model.eval(
                tile,
                diameter=params['diameter'],
                channels=params['channels'],
                flow_threshold=params['flow_threshold'],
                cellprob_threshold=params['cellprob_threshold']
            )
            
            # Convert to proper data type and update cell IDs
            masks_ij = masks_ij.astype(np.uint32)
            ncells_ij = masks_ij.max()
            
            # Shift cell IDs to avoid conflicts between tiles
            masks_ij[masks_ij != 0] += nextcell
            nextcell += ncells_ij


            # Place results into output arrays
            masks[x_from:x_to, y_from:y_to] = masks_ij
            cellprob[x_from:x_to, y_from:y_to] = flows_ij[2]  # Cell probability is in flows[2]

            if args.plot == True:
                # plot intermediate results
                print(f"Plotting {ncells_ij} cells in tile {i+1},{j+1}")
                fig = plt.figure(figsize=(12,5))
                plot.show_segmentation(fig, tile, masks_ij, flows_ij[0], channels=params['channels'])
                plt.tight_layout()
                out_name = f"cellpose_results_tile_{i+1}_{j+1}_cells{ncells_ij}.png"
                plt.savefig(out_name, dpi=300, bbox_inches='tight')
                plt.show()

            if torch.cuda.is_available():
                torch.cuda.empty_cache()
                gc.collect()
            print(f"Segmented {nextcell-ncells_ij} cells")
    print(f"Segmented {nextcell} cells")

    # Save results as compressed numpy arrays
    with gzip.open(args.masks_output, "wb") as out:
        np.save(out, masks)
    with gzip.open(args.cellprobs_output, "wb") as out:
        np.save(out, cellprob)
        
    print(f"Results saved to {args.masks_output} and {args.cellprobs_output}")