#!/bin/bash
set -euo pipefail

VCF="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/2_data/snps/Chamaea_auto_filteredQC_maf0.05_hwe0.01_prune1kb.vcf.gz"
OUTDIR="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/09_mmrr"
PLINK="/opt/homebrew/bin/plink"
DATASET="Chamaea_auto_filteredQC_maf0.05_hwe0.01_prune1kb"

BED_PREFIX="$OUTDIR/${DATASET}_mmrr"
DIST_PREFIX="$OUTDIR/${DATASET}_gendist"

mkdir -p "$OUTDIR"

"$PLINK" \
  --vcf "$VCF" \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --make-bed \
  --out "$BED_PREFIX"

"$PLINK" \
  --bfile "$BED_PREFIX" \
  --allow-extra-chr \
  --distance square 1-ibs \
  --out "$DIST_PREFIX"

cp "$DIST_PREFIX.mdist" "$OUTDIR/gendist.all.txt"
cp "$DIST_PREFIX.mdist.id" "$OUTDIR/gendist.all.ids.txt"
