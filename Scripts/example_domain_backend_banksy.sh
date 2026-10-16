#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: bash example_domain_backend_banksy.sh <input_manifest.json> <output_labels.csv>" >&2
  exit 1
fi

INPUT_MANIFEST="$1"
OUTPUT_LABELS="$2"

module load R/4.4.0-openblas-rocky8
Rscript /users/9/feng0421/Code/R/Evaluation/MainV2/Scripts/run_domain_backend_banksy.R \
  "$INPUT_MANIFEST" \
  "$OUTPUT_LABELS"
