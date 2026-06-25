#!/usr/bin/env bash
set -euo pipefail

# Post-GATK-QC filtering used to create:
# Chamaea_autoALL_filteredQC_maf0.05_miss0.25.vcf.gz
#
# Filters:
#   1. SNPs only
#   2. Biallelic variants only
#   3. Minor allele frequency (MAF) >= 0.05
#   4. No more than 25% missing genotypes (F_MISSING <= 0.25)
#
# This file is not LD-pruned and does not apply an HWE filter.
# The input VCF already contains only autosomal scaffolds ("autoALL").

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SNP_DIR="$PROJECT_DIR/2_data/snps"

INPUT_VCF="$SNP_DIR/Chamaea_autoALL_filteredQC.vcf.gz"
OUTPUT_VCF="$SNP_DIR/Chamaea_autoALL_filteredQC_maf0.05_miss0.25.vcf.gz"
LOG_FILE="$SNP_DIR/Chamaea_autoALL_filteredQC_maf0.05_miss0.25.filter.log"
EXPECTED_VARIANTS=9521177

# Prefer the working bcftools installation used in the confirmed reruns.
if [[ -x "/Users/michaeltofflemire/miniconda3/envs/vcf/bin/bcftools" ]]; then
  BCFTOOLS="/Users/michaeltofflemire/miniconda3/envs/vcf/bin/bcftools"
else
  BCFTOOLS="$(command -v bcftools || true)"
fi

if [[ -z "$BCFTOOLS" ]]; then
  echo "ERROR: bcftools was not found." >&2
  exit 1
fi

if [[ ! -s "$INPUT_VCF" ]]; then
  echo "ERROR: Input VCF is missing or empty: $INPUT_VCF" >&2
  exit 1
fi

# Avoid accidentally replacing the confirmed VCF. To intentionally rerun and
# replace it, invoke the script as: FORCE=1 bash 01b_filter_post_gatk_qc.sh
if [[ -e "$OUTPUT_VCF" && "${FORCE:-0}" != "1" ]]; then
  echo "ERROR: Output already exists: $OUTPUT_VCF" >&2
  echo "Use FORCE=1 to replace it." >&2
  exit 1
fi

TMP_VCF="${OUTPUT_VCF%.vcf.gz}.tmp.vcf.gz"
rm -f "$TMP_VCF" "$TMP_VCF.csi" "$TMP_VCF.tbi"

{
  date
  "$BCFTOOLS" --version | head -n 2
  echo "Input:  $INPUT_VCF"
  echo "Output: $OUTPUT_VCF"
  echo "Filter: SNPs; biallelic; MAF >= 0.05; F_MISSING <= 0.25"
} > "$LOG_FILE"

"$BCFTOOLS" view \
  -v snps \
  -m2 \
  -M2 \
  -i 'MAF[0]>=0.05 && F_MISSING<=0.25' \
  -Oz \
  -o "$TMP_VCF" \
  "$INPUT_VCF" \
  2>> "$LOG_FILE"

"$BCFTOOLS" index --force --tbi "$TMP_VCF" 2>> "$LOG_FILE"

VARIANT_COUNT="$("$BCFTOOLS" index -n "$TMP_VCF")"
echo "Output variants: $VARIANT_COUNT" >> "$LOG_FILE"

if [[ "$VARIANT_COUNT" -ne "$EXPECTED_VARIANTS" ]]; then
  echo "WARNING: Expected $EXPECTED_VARIANTS variants, but obtained $VARIANT_COUNT." |
    tee -a "$LOG_FILE" >&2
fi

mv "$TMP_VCF" "$OUTPUT_VCF"
mv "$TMP_VCF.tbi" "$OUTPUT_VCF.tbi"

echo "Finished: $OUTPUT_VCF"
echo "Variants retained: $VARIANT_COUNT"
echo "Log: $LOG_FILE"
