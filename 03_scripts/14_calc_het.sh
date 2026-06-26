#!/bin/bash
set -euo pipefail

VCF_FILE="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/2_data/snps/Chamaea_auto_filteredQC.vcf.gz"
RESULTS_DIR="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/06_Genome-wide_stats"
NORTH_SAMPLES="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/3_scripts/North.txt"
SOUTH_SAMPLES="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/3_scripts/South.txt"
PLINK="/opt/homebrew/bin/plink"
DATASET="Chamaea_auto_filteredQC"

BED_PREFIX="$RESULTS_DIR/${DATASET}"
HET_PREFIX="$RESULTS_DIR/${DATASET}_inbreeding_het"
SUMMARY_FILE="$RESULTS_DIR/${DATASET}_inbreeding_het_summary.tsv"
GROUP_SUMMARY_FILE="$RESULTS_DIR/${DATASET}_inbreeding_het_North_South_summary.tsv"

mkdir -p "$RESULTS_DIR"

"$PLINK" \
      --vcf "$VCF_FILE" \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --make-bed \
      --out "$BED_PREFIX"

"$PLINK" \
      --bfile "$BED_PREFIX" \
      --allow-extra-chr \
      --het \
      --out "$HET_PREFIX"

awk '
    BEGIN {
        OFS = "\t"
        print "sampleID", "observed_hom", "expected_hom", "observed_het", "expected_het", "F"
    }
    NR > 1 {
        sample = $2
        observed_hom = $3 / $5
        expected_hom = $4 / $5
        observed_het = ($5 - $3) / $5
        expected_het = ($5 - $4) / $5
        print sample, observed_hom, expected_hom, observed_het, expected_het, $6
    }
' "$HET_PREFIX.het" > "$SUMMARY_FILE"

awk '
    BEGIN {
        OFS = "\t"
        print "cluster", "n",
              "mean_observed_hom", "sd_observed_hom",
              "mean_expected_hom", "sd_expected_hom",
              "mean_observed_het", "sd_observed_het",
              "mean_expected_het", "sd_expected_het",
              "mean_F", "sd_F"
    }
    FILENAME == ARGV[1] { group[$1] = "North"; next }
    FILENAME == ARGV[2] { group[$1] = "South"; next }
    FILENAME == ARGV[3] && NR > 1 {
        cluster = group[$1]
        if (cluster == "") next

        n[cluster]++
        for (i = 2; i <= 6; i++) {
            sum[cluster, i] += $i
            sumsq[cluster, i] += $i * $i
        }
    }
    END {
        for (cluster_i = 1; cluster_i <= 2; cluster_i++) {
            cluster = cluster_i == 1 ? "North" : "South"
            printf "%s\t%d", cluster, n[cluster]
            for (i = 2; i <= 6; i++) {
                mean = sum[cluster, i] / n[cluster]
                sd = sqrt((sumsq[cluster, i] - (sum[cluster, i] * sum[cluster, i] / n[cluster])) / (n[cluster] - 1))
                printf "\t%.6f\t%.6f", mean, sd
            }
            printf "\n"
        }
    }
' "$NORTH_SAMPLES" "$SOUTH_SAMPLES" "$SUMMARY_FILE" > "$GROUP_SUMMARY_FILE"
