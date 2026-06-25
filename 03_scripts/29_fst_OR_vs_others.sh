#!/bin/bash

set -euo pipefail

VCF_FILE="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/data/Chamaea_autoALL_filteredQC_maf0.05_miss0.25.vcf.gz"

POP_DIR="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts"
OREGON_POP_FILE="$POP_DIR/OC.txt"

RESULTS_DIR="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Selective_sweeps"
mkdir -p "$RESULTS_DIR"

for POP_FILE in "$POP_DIR"/*.txt; do
    if [[ "$POP_FILE" == "$OREGON_POP_FILE" ]]; then
        continue
    fi

    POP_NAME=$(basename "$POP_FILE" .txt)
    FST_OUTPUT_PREFIX="$RESULTS_DIR/OC_vs_${POP_NAME}_fst"

    echo "Running Fst: OC vs ${POP_NAME}"

    vcftools --gzvcf "$VCF_FILE" \
             --weir-fst-pop "$OREGON_POP_FILE" \
             --weir-fst-pop "$POP_FILE" \
             --fst-window-size 50000 \
             --fst-window-step 10000 \
             --out "$FST_OUTPUT_PREFIX"

    echo "Finished: OC vs ${POP_NAME}"
done

echo "All Oregon pairwise Fst analyses complete."
echo "Results saved in: $RESULTS_DIR"