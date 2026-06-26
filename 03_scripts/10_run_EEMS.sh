#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

OUTDIR="$PROJECT_DIR/4_analyses/05_eems"
PREFIX="Chamaea_auto_filteredQC_eems"
INI="$OUTDIR/$PREFIX.ini"
EEMS="/Users/michaeltofflemire/softs/eems_clean/runeems_snps/src/runeems_snps"
SEED=123
LOG="$OUTDIR/${PREFIX}_runeems_seed${SEED}.log"

"$EEMS" --params "$INI" --seed "$SEED" > "$LOG" 2>&1
