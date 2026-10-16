import argparse
from pathlib import Path

import pandas as pd


def read_values_from_file(file_path):
    """Read non-empty values from a text file (one per line)."""
    with open(file_path, "r", encoding="utf-8") as f:
        values = [line.strip() for line in f if line.strip()]
    return values


def default_generated_targets_path(input_csv):
    path = Path(input_csv)
    name = path.name
    for suffix in [".csv.gz", ".csv"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return str(path.with_name(f"{name}_target_to_filter.txt"))


def derive_targets_from_prefixes(df, target_column, prefixes):
    series = df[target_column].dropna().astype(str)
    matched = sorted({value for value in series if any(value.startswith(prefix) for prefix in prefixes)})
    return matched


def write_values_to_file(file_path, values):
    with open(file_path, "w", encoding="utf-8") as f:
        for value in values:
            f.write(f"{value}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Filter out rows where target column values are in the specified list or match configured prefixes."
    )
    parser.add_argument("input_csv", help="Input CSV file path")
    parser.add_argument("output_csv", help="Output CSV file path")

    parser.add_argument(
        "--targets",
        "-t",
        nargs="+",
        help="List of target values to filter out",
    )
    parser.add_argument(
        "--targets-file",
        "-f",
        help="File containing target values to filter out (one per line)",
    )
    parser.add_argument(
        "--prefixes",
        nargs="+",
        help="Filter out all values whose target column starts with any of these prefixes",
    )
    parser.add_argument(
        "--prefixes-file",
        help="File containing prefixes to match (one per line)",
    )
    parser.add_argument(
        "--generated-targets-file",
        help="Optional path for writing the targets discovered via prefix matching; defaults to <input>_target_to_filter.txt",
    )
    parser.add_argument(
        "--target-column",
        default="target",
        help="Name of the target column (default: target)",
    )

    args = parser.parse_args()

    explicit_target_flags = sum(bool(option) for option in [args.targets, args.targets_file])
    prefix_flags = sum(bool(option) for option in [args.prefixes, args.prefixes_file])
    if explicit_target_flags > 1:
        parser.error("Use only one of --targets or --targets-file.")
    if prefix_flags > 1:
        parser.error("Use only one of --prefixes or --prefixes-file.")
    if explicit_target_flags == 0 and prefix_flags == 0:
        parser.error("Provide one of --targets, --targets-file, --prefixes, or --prefixes-file.")

    print(f"Reading {args.input_csv}...")
    df = pd.read_csv(args.input_csv)

    targets_to_filter = []

    if args.targets_file:
        targets_to_filter = read_values_from_file(args.targets_file)
        print(f"Read {len(targets_to_filter)} targets from {args.targets_file}")
        print(targets_to_filter[:10])
    elif args.targets:
        targets_to_filter = args.targets

    prefixes = []
    if args.prefixes_file:
        prefixes = read_values_from_file(args.prefixes_file)
        print(f"Read {len(prefixes)} prefixes from {args.prefixes_file}")
        print(prefixes[:10])
    elif args.prefixes:
        prefixes = args.prefixes

    if prefixes:
        prefix_targets = derive_targets_from_prefixes(df, args.target_column, prefixes)
        generated_targets_file = args.generated_targets_file or default_generated_targets_path(args.input_csv)
        write_values_to_file(generated_targets_file, prefix_targets)
        print(f"Generated {len(prefix_targets)} targets from prefixes and saved to: {generated_targets_file}")
        print(prefix_targets[:10])
        targets_to_filter = sorted(set(targets_to_filter).union(prefix_targets))

    if not targets_to_filter:
        print("No targets matched; output will be identical to input.")

    df_filtered = df[~df[args.target_column].isin(targets_to_filter)]

    print(f"Original rows: {len(df)}, Filtered rows: {len(df_filtered)}, Removed: {len(df) - len(df_filtered)}")
    print(f"Targets to filter: {len(targets_to_filter)} items")

    df_filtered.to_csv(args.output_csv, index=False)
    print(f"Filtered data saved to: {args.output_csv}")


if __name__ == "__main__":
    main()


# Example usage:
# cd /users/9/feng0421/Code/BIDCell/sc_ref
# python clean_trans.py /users/9/feng0421/BIDCell/data/MerFish/HumanLungCancerPatient1_detected_transcripts.csv /users/9/feng0421/BIDCell/data/MerFish/HumanLungCancerPatient1_detected_transcripts_filtered.csv -f /users/9/feng0421/Code/BIDCell/sc_ref/gene_to_filter.txt --target-column gene
# python clean_trans.py /users/9/feng0421/BIDCell/data/MerFish/roi_transcripts.csv /users/9/feng0421/BIDCell/data/MerFish/roi_detected_transcripts_filtered.csv -f /users/9/feng0421/Code/BIDCell/sc_ref/gene_to_filter.txt --target-column gene
# python clean_trans.py /users/9/feng0421/BIDCell/data/CosMx/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_tx_file.csv /scratch.global/feng0421/data/CosMx/Lung5_Rep1_tx_file.csv -f /users/9/feng0421/Code/BIDCell/sc_ref/target_to_filter.txt --target-column target
# python clean_trans.py /scratch.global/feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts.csv.gz /scratch.global/feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts_filtered.csv --target-column feature_name --prefixes NegControlProbe_ NegControlCodeword_ antisense_ BLANK_ Blank- NegPrb
