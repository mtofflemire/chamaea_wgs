############################################################
## RDA/LFMM overlapping genes
## Uses current final files:
## RDA:  results/10_rda/RDA_outlier_GENES.csv
## LFMM: results/11_lfmm/LFMM_outlier_GENES.csv
############################################################

library(dplyr)

setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea")

rda_gene_hits <- read.csv(
  "results/10_rda/RDA_outlier_GENES.csv",
  stringsAsFactors = FALSE
)

lfmm_gene_hits <- read.csv(
  "results/11_lfmm/LFMM_outlier_GENES.csv",
  stringsAsFactors = FALSE
)

############################################################
## Make unique gene-level tables
############################################################

rda_genes_unique <- rda_gene_hits %>%
  dplyr::filter(!is.na(gene_id), gene_id != "") %>%
  dplyr::distinct(gene_id, gene_name, gene_biotype, description) %>%
  dplyr::arrange(gene_name)

lfmm_genes_unique <- lfmm_gene_hits %>%
  dplyr::filter(!is.na(gene_id), gene_id != "") %>%
  dplyr::distinct(gene_id, gene_name, gene_biotype, description) %>%
  dplyr::arrange(gene_name)

############################################################
## Find overlapping genes
############################################################

rda_lfmm_overlapping_genes <- rda_genes_unique %>%
  dplyr::inner_join(
    lfmm_genes_unique,
    by = c("gene_id", "gene_name", "gene_biotype", "description")
  ) %>%
  dplyr::arrange(gene_name)

############################################################
## Add SNP positions from each method
############################################################

rda_gene_positions <- rda_gene_hits %>%
  dplyr::filter(gene_id %in% rda_lfmm_overlapping_genes$gene_id) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    rda_snps = paste(sort(unique(SNP)), collapse = "; "),
    rda_snp_positions = paste(sort(unique(position)), collapse = "; "),
    n_rda_snps = dplyr::n_distinct(SNP),
    .groups = "drop"
  )

lfmm_gene_positions <- lfmm_gene_hits %>%
  dplyr::filter(gene_id %in% rda_lfmm_overlapping_genes$gene_id) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    lfmm_snps = paste(sort(unique(snp)), collapse = "; "),
    lfmm_snp_positions = paste(sort(unique(position)), collapse = "; "),
    lfmm_variables = paste(sort(unique(var)), collapse = "; "),
    n_lfmm_snps = dplyr::n_distinct(snp),
    .groups = "drop"
  )

rda_lfmm_overlapping_genes_full <- rda_lfmm_overlapping_genes %>%
  dplyr::left_join(rda_gene_positions, by = "gene_id") %>%
  dplyr::left_join(lfmm_gene_positions, by = "gene_id") %>%
  dplyr::arrange(gene_name)

############################################################
## Save outputs
############################################################



write.csv(
  rda_lfmm_overlapping_genes,
  "results/12_GO/RDA_LFMM_overlapping_GENES.csv",
  row.names = FALSE
)

write.csv(
  rda_lfmm_overlapping_genes_full,
  "results/12_GO/RDA_LFMM_overlapping_GENES_full.csv",
  row.names = FALSE
)

writeLines(
  rda_lfmm_overlapping_genes$gene_name,
  "results/12_GO/RDA_LFMM_overlapping_GENES.txt"
)

############################################################
## Print summary
############################################################

cat("Unique RDA genes:", nrow(rda_genes_unique), "\n")
cat("Unique LFMM genes:", nrow(lfmm_genes_unique), "\n")
cat("Overlapping RDA/LFMM genes:", nrow(rda_lfmm_overlapping_genes), "\n")


















# Use this after downloading the PANTHER GO results table

# Install once if needed:
# install.packages("BiocManager")
BiocManager::install(c("GO.db", "org.Gg.eg.db"))

library(readr)
library(stringr)
library(dplyr)
library(GO.db)
library(AnnotationDbi)
library(org.Gg.eg.db)

analysis_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/12_go/analysis.txt"

# CHANGE THIS to your original uploaded gene list
# It should be one gene ID per line
gene_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/12_go/upload_1.txt"

# CHANGE THIS depending on your gene IDs:
# Common chicken options: "SYMBOL", "ENSEMBL", "ENTREZID"
gene_keytype <- "SYMBOL"

# Read PANTHER output
lines <- readLines(analysis_file)

header_i <- grep("^GO biological process complete", lines)

go_raw <- readr::read_tsv(
  paste(lines[header_i:length(lines)], collapse = "\n"),
  show_col_types = FALSE
)

colnames(go_raw) <- c(
  "GO_term_full",
  "Reference_count",
  "Input_count",
  "Expected",
  "Direction",
  "Fold_enrichment",
  "Raw_p_value",
  "FDR"
)

# Clean GO enrichment table
go_table <- go_raw %>%
  dplyr::mutate(
    GO_ID = stringr::str_extract(GO_term_full, "GO:\\d+|UNCLASSIFIED"),
    GO_term = stringr::str_trim(
      stringr::str_remove(GO_term_full, "\\s*\\((GO:\\d+|UNCLASSIFIED)\\)$")
    )
  ) %>%
  dplyr::select(
    GO_term,
    GO_ID,
    Reference_count,
    Input_count,
    Expected,
    Direction,
    Fold_enrichment,
    Raw_p_value,
    FDR
  )

# Get official GO term names and definitions
go_definitions <- AnnotationDbi::select(
  GO.db::GO.db,
  keys = go_table$GO_ID[go_table$GO_ID != "UNCLASSIFIED"],
  columns = c("TERM", "DEFINITION", "ONTOLOGY"),
  keytype = "GOID"
)

go_table_annotated <- go_table %>%
  dplyr::left_join(go_definitions, by = c("GO_ID" = "GOID")) %>%
  dplyr::mutate(
    TERM = dplyr::if_else(is.na(TERM), GO_term, TERM),
    DEFINITION = dplyr::if_else(
      GO_ID == "UNCLASSIFIED",
      "Unclassified genes without GO assignment",
      DEFINITION
    )
  )

# Read your uploaded/input genes
input_genes <- readr::read_lines(gene_file)
input_genes <- stringr::str_trim(input_genes)
input_genes <- input_genes[input_genes != ""]

# Map input genes to GO terms
gene_go_map <- AnnotationDbi::select(
  org.Gg.eg.db::org.Gg.eg.db,
  keys = input_genes,
  keytype = gene_keytype,
  columns = c(gene_keytype, "GOALL", "ONTOLOGYALL")
)

gene_go_map <- gene_go_map %>%
  dplyr::filter(
    !is.na(GOALL),
    ONTOLOGYALL == "BP"
  ) %>%
  dplyr::rename(
    Input_gene = !!gene_keytype,
    GO_ID = GOALL
  )

# Collapse genes for each GO term
genes_by_go <- gene_go_map %>%
  dplyr::group_by(GO_ID) %>%
  dplyr::summarise(
    Associated_input_genes = paste(sort(unique(Input_gene)), collapse = "; "),
    Associated_input_gene_count = dplyr::n_distinct(Input_gene),
    .groups = "drop"
  )

# Add genes to annotated GO table
go_table_with_genes <- go_table_annotated %>%
  dplyr::left_join(genes_by_go, by = "GO_ID") %>%
  dplyr::mutate(
    Associated_input_genes = dplyr::if_else(
      is.na(Associated_input_genes),
      "",
      Associated_input_genes
    ),
    Associated_input_gene_count = dplyr::if_else(
      is.na(Associated_input_gene_count),
      0L,
      Associated_input_gene_count
    )
  )

# Save output
output_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/12_go/go_terms_annotated_with_genes.tsv"

readr::write_tsv(go_table_with_genes, output_file)

# View top 20 enriched GO terms
go_table_with_genes %>%
  dplyr::arrange(FDR) %>%
  head(20)

View(go_table_with_genes)