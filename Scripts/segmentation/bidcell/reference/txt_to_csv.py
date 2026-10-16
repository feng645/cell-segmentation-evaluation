import pandas as pd
import io, re

in_txt  = "/scratch.global/feng0421/data/ref/NSCLC_GSE148071_expression_Celltype_majorlineage.txt"      # input txt file; update the atlas name below accordingly
out_csv = "/scratch.global/feng0421/data/ref/NSCLC_GSE148071_expression_Celltype_majorlineage.csv"

# 1) Detect the delimiter (tab or multiple spaces)
with open(in_txt, "r", encoding="utf-8") as f:
    head = "".join([next(f) for _ in range(5)])
sep = "\t" if "\t" in head else r"\s+"

# 2) Read: first column = gene names (index), other columns = cell types
df = pd.read_csv(in_txt, sep=sep, engine="python", header=0, index_col=0)


# 4) Average duplicated gene symbols
#    (gene names are in the index, so aggregate over axis=0)
df = df.groupby(level=0).mean()

# 5) Transpose: cell types -> rows, genes -> columns
df = df.T


ct_map = {
    "AT2": 0,
    "AT1": 1,
    "Smooth muscle": 2,
    "CD8 T cells": 3,
    "CD4 T cells": 4,
    "T cells proliferating": 5,
    "NK cells": 6,
    "B cells": 7,
    "Plasma cells": 8,
    "Classical monocytes": 9,
    "Alveolar macrophages": 10,
    "Monocyte-derived Mph": 11,
    "Fibroblasts": 12,
    "Endothelial": 13,
    "Malignant": 14,
}

# Keep only the cell types in ct_map
ct_keys = list(ct_map.keys())
df = df[df.index.isin(ct_keys)]

df.insert(len(df.columns), "ct_idx", [ct_map[ct] for ct in df.index])
df.insert(len(df.columns), "cell_type", df.index)
df.insert(len(df.columns), "atlas", "NSCLC_GSE148071")

# 6) Save, keeping the row index (cell-type names)
df.to_csv(out_csv, index=True)
print(f"Done: {out_csv}\nshape={df.shape}  (rows=cells, cols=genes)")