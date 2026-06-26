#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

VCF="$PROJECT_DIR/2_data/snps/Chamaea_auto_filteredQC.vcf.gz"
OUTDIR="$PROJECT_DIR/4_analyses/05_eems"
PREFIX="Chamaea_auto_filteredQC_eems"
BED_PREFIX="$OUTDIR/$PREFIX"
DATAPATH="$OUTDIR/$PREFIX"
MCMC_DIR="$OUTDIR/chamaea_eems_run1"

METADATA="$PROJECT_DIR/1_Meta/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv"
OUTER_TEMPLATE="$SCRIPT_DIR/eems.outer"
MAKE_COORD_SCRIPT="$SCRIPT_DIR/08_make_coord_mat.R"
INI="$OUTDIR/$PREFIX.ini"

mkdir -p "$OUTDIR" "$MCMC_DIR"

plink \
  --vcf "$VCF" \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --make-bed \
  --threads 8 \
  --out "$BED_PREFIX"

plink \
  --bfile "$BED_PREFIX" \
  --allow-extra-chr \
  --distance square 1-ibs \
  --threads 8 \
  --out "$BED_PREFIX"

cp "$BED_PREFIX.mdist" "$DATAPATH.diffs"
cp "$OUTER_TEMPLATE" "$DATAPATH.outer"
Rscript "$MAKE_COORD_SCRIPT"

NSITES=$(wc -l < "$BED_PREFIX.bim" | tr -d ' ')
NINDIV=$(wc -l < "$BED_PREFIX.fam" | tr -d ' ')

cat > "$INI" <<EOF
datapath = $DATAPATH
mcmcpath = $MCMC_DIR
nIndiv = $NINDIV
nSites = $NSITES
nDemes = 400
diploid = true
numMCMCIter = 2000000
numBurnIter = 1000000
numThinIter = 9999
EOF

cat > "$OUTDIR/${PREFIX}_input_summary.txt" <<EOF
VCF: $VCF
PLINK prefix: $BED_PREFIX
EEMS datapath prefix: $DATAPATH
EEMS mcmcpath: $MCMC_DIR
EEMS ini: $INI
nIndiv: $NINDIV
nSites: $NSITES
outer template: $OUTER_TEMPLATE
metadata: $METADATA
EOF
