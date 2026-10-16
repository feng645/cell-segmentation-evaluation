# since BIDCell requires gene in reference files must match those in transcript files, 
# so this file will extract gene list from transcript files

import pandas as pd

# Input CSV file path and output TXT file path
in_csv = "/scratch.global/feng0421/data/MERFISH/dataset/HumanLungCancerPatient1_detected_transcripts.csv"
out_txt = "/scratch.global/feng0421/data/MERFISH/dataset/Merfish_gene.txt"

# Read the CSV (reading only the target column is faster)
df = pd.read_csv(in_csv, usecols=["gene"])

# Unique values
unique_targets = df["gene"].dropna().unique()

# Save to a text file
with open(out_txt, "w") as f:
    for val in sorted(unique_targets):
        f.write(str(val) + "\n")

print(f"Unique values from 'gene' column saved to {out_txt}, total={len(unique_targets)}")
