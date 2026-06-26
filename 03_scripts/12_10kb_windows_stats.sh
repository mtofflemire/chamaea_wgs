#!/bin/bash
set -euo pipefail

VCF_FILE="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/2_data/snps/Chamaea_auto_filteredQC.vcf.gz"
POP_NORTH_FILE="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/3_scripts/North.txt"
POP_SOUTH_FILE="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/3_scripts/South.txt"
RESULTS_DIR="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/06_Genome-wide_stats"
VCFTOOLS="/Users/michaeltofflemire/miniconda3/envs/vcf/bin/vcftools"
DATASET="Chamaea_auto_filteredQC"

mkdir -p "$RESULTS_DIR"

"$VCFTOOLS" --gzvcf "$VCF_FILE" \
         --weir-fst-pop "$POP_NORTH_FILE" \
         --weir-fst-pop "$POP_SOUTH_FILE" \
         --fst-window-size 10000 \
         --fst-window-step 10000 \
         --out "$RESULTS_DIR/${DATASET}_North_vs_South_10kb"

"$VCFTOOLS" --gzvcf "$VCF_FILE" \
         --keep "$POP_NORTH_FILE" \
         --window-pi 10000 \
         --out "$RESULTS_DIR/${DATASET}_North_10kb"

"$VCFTOOLS" --gzvcf "$VCF_FILE" \
         --keep "$POP_SOUTH_FILE" \
         --window-pi 10000 \
         --out "$RESULTS_DIR/${DATASET}_South_10kb"

"$VCFTOOLS" --gzvcf "$VCF_FILE" \
         --keep "$POP_NORTH_FILE" \
         --TajimaD 10000 \
         --out "$RESULTS_DIR/${DATASET}_North_10kb"

"$VCFTOOLS" --gzvcf "$VCF_FILE" \
         --keep "$POP_SOUTH_FILE" \
         --TajimaD 10000 \
         --out "$RESULTS_DIR/${DATASET}_South_10kb"
