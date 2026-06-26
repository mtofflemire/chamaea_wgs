#!/usr/bin/env Rscript

sample_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/3_scripts/samples.txt"
metadata_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/1_Meta/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv"
out_coord <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/05_eems/Chamaea_auto_filteredQC_eems.coord"

samples <- read.table(sample_file, stringsAsFactors = FALSE)[, 1]
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

coords <- meta[match(samples, meta$sampleID), c("long", "lat")]

write.table(
    coords,
    file = out_coord,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)
