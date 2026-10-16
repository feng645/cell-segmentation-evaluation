from PIL import Image
import numpy as np
import json
import cv2

def tif_to_geojson(tif_path, output_path, simplify_tolerance=1.0):
    """
    Convert a TIF cell-label image to GeoJSON.
    
    Parameters:
    - tif_path: path to the TIF file
    - output_path: path to the output JSON file  
    - simplify_tolerance: contour simplification tolerance (fewer vertices)
    """
    
    # Read the TIF image
    print("Reading TIF file...")
    img = Image.open(tif_path)
    label_image = np.array(img)
    
    print(f"Image size: {label_image.shape}")
    print(f"Pixel value range: {label_image.min()} - {label_image.max()}")
    
    # Collect all non-background cell labels
    unique_labels = np.unique(label_image)
    unique_labels = unique_labels[unique_labels != 0]  # exclude background (0)
    
    print(f"Detected {len(unique_labels)} cells")
    
    features = []
    
    # Process each cell
    for i, cell_id in enumerate(unique_labels):
        if i % 100 == 0:  # report progress every 100 cells
            print(f"Progress: {i}/{len(unique_labels)}")
            
        # Binary mask of the current cell
        mask = (label_image == cell_id).astype(np.uint8)
        
        mask = mask * 255

        # Find contours and their hierarchy with OpenCV, keeping the ring structure
        contours_info = cv2.findContours(mask, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
        if len(contours_info) == 3:
            _img, contours, hierarchy = contours_info
        else:
            contours, hierarchy = contours_info

        if hierarchy is None or len(contours) == 0:
            continue

        hierarchy = hierarchy[0]

        def contour_to_ring(contour_points):
            coords_array = contour_points.reshape(-1, 1, 2).astype(np.float32)
            ring = coords_array.copy()
            if simplify_tolerance > 0 and len(coords_array) > 2:
                arc_length = cv2.arcLength(coords_array, True)
                epsilon = simplify_tolerance if simplify_tolerance >= 1 else simplify_tolerance * arc_length
                simplified = cv2.approxPolyDP(coords_array, epsilon, True)
                if len(simplified) >= 3:
                    ring = simplified
            coords = [[float(p[0][0]), float(p[0][1])] for p in ring]
            if coords and coords[0] != coords[-1]:
                coords.append(coords[0])
            return coords if len(coords) >= 4 else []

        polygons = []

        for idx, contour in enumerate(contours):
            parent_idx = hierarchy[idx][3]
            if parent_idx != -1:
                continue  # outer contours only; holes are attached below

            outer_ring = contour_to_ring(contour)
            if not outer_ring:
                continue

            holes = []
            child_idx = hierarchy[idx][2]
            while child_idx != -1:
                hole_ring = contour_to_ring(contours[child_idx])
                if hole_ring:
                    holes.append(hole_ring)
                child_idx = hierarchy[child_idx][0]

            polygons.append([outer_ring] + holes)

        if not polygons:
            continue

        if len(polygons) == 1:
            geometry = {
                "type": "Polygon",
                "coordinates": polygons[0]
            }
        else:
            geometry = {
                "type": "MultiPolygon",
                "coordinates": polygons
            }

        feature = {
            "geometry": geometry,
            "id": str(int(cell_id)),
            "type": "Feature"
        }
        features.append(feature)
    
    # Build the GeoJSON structure
    geojson = {
        "features": features,
        "type": "FeatureCollection"
    }
    
    # Save the JSON file
    print(f"Saved to: {output_path}")
    with open(output_path, 'w') as f:
        json.dump(geojson, f, indent=2)
    
    print(f"Done: {len(features)} cell polygons")
    return geojson

# Example usage
if __name__ == "__main__":
    # Input TIF path
    tif_file = "/projects/standard/zhangl4/shared/Data_feng0421/Xenium_BIDCell/model_outputs/2026_04_04_00_11_19/test_output/epoch_1_step_4000_connected.tif"
    
    # Output JSON path
    output_file = "/projects/standard/zhangl4/shared/Data_feng0421/Xenium_BIDCell/model_outputs/2026_04_04_00_11_19/test_output/cells_output.json"
    
    # Run the conversion
    # simplify_tolerance: 0 disables simplification; 1-2 is a reasonable amount
    result = tif_to_geojson(tif_file, output_file, simplify_tolerance=1.0)
    
    # Print a few statistics
    print(f"\nThe JSON contains {len(result['features'])} cells")
    if len(result['features']) > 0:
        first_feature = result['features'][0]
        coord_count = len(first_feature['geometry']['coordinates'][0])
        print(f"Vertices in the first cell: {coord_count}")
        print(f"First cell ID: {first_feature['id']}")
