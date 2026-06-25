#!/usr/bin/env bash
set -euo pipefail

# Windowed Weir and Cockerham Fst between the northern and southern
# ADMIXTURE-defined genetic groups.
#
# Input: unpruned autosomal VCF (9,521,177 SNPs in the confirmed analysis)
# Windows: 50 kb, advanced in 10-kb steps

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

VCF="$PROJECT_DIR/2_data/snps/Chamaea_autoALL_filteredQC_maf0.05_miss0.25.vcf.gz"
NORTH_SAMPLES="$SCRIPT_DIR/north_genetic_group_samples.txt"
SOUTH_SAMPLES="$SCRIPT_DIR/south_genetic_group_samples.txt"
OUTPUT_DIR="$PROJECT_DIR/4_results/Selective_sweeps/North_vs_South_Fst"
OUTPUT_PREFIX="$OUTPUT_DIR/North_vs_South_50kb_step10kb"

if [[ -x "/Users/michaeltofflemire/miniconda3/envs/vcf/bin/vcftools" ]]; then
  VCFTOOLS="/Users/michaeltofflemire/miniconda3/envs/vcf/bin/vcftools"
else
  VCFTOOLS="$(command -v vcftools || true)"
fi

if [[ -z "$VCFTOOLS" ]]; then
  echo "ERROR: vcftools was not found." >&2
  exit 1
fi

for input_file in "$VCF" "$NORTH_SAMPLES" "$SOUTH_SAMPLES"; do
  if [[ ! -s "$input_file" ]]; then
    echo "ERROR: Missing or empty input file: $input_file" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

echo "Northern samples: $(wc -l < "$NORTH_SAMPLES" | tr -d ' ')"
echo "Southern samples: $(wc -l < "$SOUTH_SAMPLES" | tr -d ' ')"
echo "Running 50-kb windowed Fst with a 10-kb step..."

"$VCFTOOLS" \
  --gzvcf "$VCF" \
  --weir-fst-pop "$NORTH_SAMPLES" \
  --weir-fst-pop "$SOUTH_SAMPLES" \
  --fst-window-size 50000 \
  --fst-window-step 10000 \
  --out "$OUTPUT_PREFIX" \
  > "$OUTPUT_DIR/North_vs_South_Fst.stdout.log" \
  2> "$OUTPUT_DIR/North_vs_South_Fst.stderr.log"

if [[ ! -s "$OUTPUT_PREFIX.windowed.weir.fst" ]]; then
  echo "ERROR: Fst output was not created." >&2
  exit 1
fi

echo "Finished."
echo "Output: $OUTPUT_PREFIX.windowed.weir.fst"
