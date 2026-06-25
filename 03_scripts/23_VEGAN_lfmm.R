############################################################
## LFMM analysis
## Uses same SNPs and environmental PCs as RDA
############################################################

library(algatr)
library(vcfR)
library(dplyr)
library(ggplot2)

setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea")

dir.create("results/11_lfmm", recursive = TRUE, showWarnings = FALSE)

############################################################
## Load sample order
############################################################

sample_ids <- readLines("scripts/samples.txt")

############################################################
## Load genotype matrix if needed
############################################################

if (!exists("genotypes")) {
  vcf <- read.vcfR("results/01_vcfs/GATK_QC/AUTOSOMES/Chamaea_autoALL_filteredQC_maf0.05_miss0.25_prune1kb.vcf.gz")
  
  gt_raw <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  
  gt_raw[gt_raw %in% c(".", "./.", ".|.")] <- NA
  gt_raw[gt_raw %in% c("0/0", "0|0")] <- 0
  gt_raw[gt_raw %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  gt_raw[gt_raw %in% c("1/1", "1|1")] <- 2
  
  genotypes <- gt_raw
  storage.mode(genotypes) <- "numeric"
  
  genotypes <- t(genotypes)
  genotypes <- genotypes[sample_ids, , drop = FALSE]
  
  genotypes <- apply(genotypes, 2, function(x) {
    replace(x, is.na(x), as.numeric(names(which.max(table(x)))))
  })
  
  genotypes <- as.matrix(genotypes)
  rownames(genotypes) <- sample_ids
}

cat("Genotype matrix dimensions:", dim(genotypes), "\n")

############################################################
## Environmental predictors
############################################################

env <- as.data.frame(env_pc_fixed)
rownames(env) <- sample_ids

cat("Environmental matrix dimensions:", dim(env), "\n")
print(colnames(env))

############################################################
## Run LFMM
############################################################

K <- 2

ridge_results <- lfmm_run(
  gen = genotypes,
  env = env,
  K = K,
  lfmm_method = "ridge",
  p_adj = "fdr",
  sig = 0.05
)



saveRDS(
  ridge_results,
  "results/12_lfmm/lfmm_ridge_results_K2.rds"
)



############################################################
## Save full LFMM result tables
############################################################

write.csv(
  ridge_results$df,
  "results/11_lfmm/LFMM_results.csv",
  row.names = FALSE
)




############################################################
## Save significant SNPs
############################################################



ridge_snps <- ridge_results$lfmm_snps


write.csv(
  ridge_snps,
  "results/11_lfmm/lfmm_outlier_SNPs.csv",
  row.names = FALSE
)



cat("Ridge significant SNPs:", nrow(ridge_snps), "\n")

ridge_snps <- ridge_results$lfmm_snps


cat("Ridge significant SNP-variable hits:", nrow(ridge_snps), "\n")




############################################################
## Add scaffold and position info
############################################################

ridge_snps <- ridge_snps %>%
  dplyr::mutate(
    scaffold = sub("_[^_]+$", "", snp),
    position = as.numeric(sub("^.*_", "", snp))
  )



write.csv(
  ridge_snps,
  "results/11_lfmm/lfmm_outlier_SNPs.csv",
  row.names = FALSE
)



############################################################
## Count SNPs by environmental PC
############################################################

ridge_snps_by_variable <- ridge_snps %>%
  dplyr::group_by(var) %>%
  dplyr::summarise(num_snps = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(desc(num_snps))






print(ridge_snps_by_variable)




























############################################################
## Manhattan plots for LFMM ridge, faceted by environmental PC
############################################################

library(dplyr)
library(ggplot2)

lfmm_df <- ridge_results$df

lfmm_df <- lfmm_df %>%
  dplyr::mutate(
    scaffold = sub("_[^_]+$", "", snp),
    snp_position = as.numeric(sub("^.*_", "", snp))
  )

scaffold_offsets <- lfmm_df %>%
  dplyr::distinct(scaffold, snp_position) %>%
  dplyr::group_by(scaffold) %>%
  dplyr::summarise(scaffold_length = max(snp_position, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(scaffold) %>%
  dplyr::mutate(
    scaffold_offset = cumsum(dplyr::lag(scaffold_length, default = 0)),
    scaffold_index = dplyr::row_number()
  )

lfmm_df <- lfmm_df %>%
  dplyr::left_join(scaffold_offsets, by = "scaffold") %>%
  dplyr::mutate(
    cumulative_position = snp_position + scaffold_offset,
    scaffold_color = ifelse(scaffold_index %% 2 == 0, "even", "odd"),
    significant = adjusted.pvalue < 0.1
  )








############################################################
## Manhattan plots for LFMM ridge
## One plot per environmental PC, no facets
############################################################

library(dplyr)
library(ggplot2)

lfmm_variables <- sort(unique(lfmm_df$var))


manhattan_plot_ridge <- ggplot(
  lfmm_df,
  aes(x = cumulative_position, y = -log10(adjusted.pvalue), color = scaffold_color)
) +
  geom_point(alpha = 0.7, size = 0.35) +
  geom_hline(
    yintercept = -log10(0.1),
    color = "red",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c("odd" = "orange2", "even" = "navyblue")) +
  facet_wrap(~ var, ncol = 2, scales = "free_x") +
  theme_bw() +
  labs(
    x = NULL,
    y = "-log10(FDR-adjusted p-value)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(manhattan_plot_ridge)







install.packages("ggh4x")
library(ggh4x)

manhattan_plot_ridge <- ggplot(
  lfmm_df,
  aes(x = cumulative_position, y = -log10(adjusted.pvalue), color = scaffold_color)
) +
  geom_point(alpha = 0.7, size = 0.35) +
  geom_hline(
    yintercept = -log10(0.1),
    color = "red",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c("odd" = "orange2", "even" = "navyblue")) +
  facet_wrap2(~ var, ncol = 2, scales = "free", axes = "all") +
  theme_bw() +
  labs(
    x = NULL,
    y = "-log10(FDR-adjusted p-value)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.text.y = element_text(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )



ggsave(
  "results/11_lfmm/lfmm_ridge_manhattan_by_variable_K2_FDR0.1.pdf",
  manhattan_plot_ridge,
  width = 10,
  height = 6
)

ggsave(
  "results/11_lfmm/lfmm_ridge_manhattan_by_variable_K2_FDR0.1.png",
  manhattan_plot_ridge,
  width = 10,
  height = 6,
  dpi = 600
)



print(manhattan_plot_ridge)

























############################################################
## LFMM ridge SNPs: inside-gene annotation only
############################################################

library(jsonlite)
library(dplyr)
library(readr)
library(stringr)

setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea")

# use your ridge significant SNP table
lfmm_outliers <- ridge_snps

# if needed, reload from saved file instead:
# lfmm_outliers <- read.csv("results/12_lfmm/lfmm_ridge_significant_snps_K2_FDR0.1_with_positions.csv")

lfmm_outliers$scaffold <- sub("_[^_]+$", "", lfmm_outliers$snp)
lfmm_outliers$position <- as.numeric(sub("^.*_", "", lfmm_outliers$snp))

seq_report <- jsonlite::stream_in(
  file("data/reference_annotation/ncbi_dataset/data/GCF_029207755.1/sequence_report.jsonl"),
  verbose = FALSE
)

name_map <- seq_report %>%
  dplyr::select(
    scaffold = genbankAccession,
    annotated_scaffold = refseqAccession
  )

lfmm_outliers_mapped <- lfmm_outliers %>%
  dplyr::left_join(name_map, by = "scaffold")

gff <- readr::read_tsv(
  "data/reference_annotation/ncbi_dataset/data/GCF_029207755.1/genomic.gff",
  comment = "#",
  col_names = c(
    "seqid", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes"
  ),
  col_types = readr::cols(
    seqid = readr::col_character(),
    source = readr::col_character(),
    type = readr::col_character(),
    start = readr::col_double(),
    end = readr::col_double(),
    score = readr::col_character(),
    strand = readr::col_character(),
    phase = readr::col_character(),
    attributes = readr::col_character()
  )
)

genes <- gff %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(
    gene_id = stringr::str_match(attributes, "ID=([^;]+)")[, 2],
    gene_name = stringr::str_match(attributes, "Name=([^;]+)")[, 2],
    gene_biotype = stringr::str_match(attributes, "gene_biotype=([^;]+)")[, 2],
    description = stringr::str_match(attributes, "description=([^;]+)")[, 2]
  )

inside_hits <- list()

for (i in seq_len(nrow(lfmm_outliers_mapped))) {
  snp <- lfmm_outliers_mapped[i, ]
  
  hits <- genes %>%
    dplyr::filter(
      seqid == snp$annotated_scaffold,
      start <= snp$position,
      end >= snp$position
    )
  
  if (nrow(hits) > 0) {
    hits$snp <- snp$snp
    hits$var <- snp$var
    hits$scaffold <- snp$scaffold
    hits$position <- snp$position
    hits$annotated_scaffold <- snp$annotated_scaffold
    inside_hits[[length(inside_hits) + 1]] <- hits
  }
}

lfmm_inside_genes_clean <- dplyr::bind_rows(inside_hits)

lfmm_snp_gene_table <- lfmm_inside_genes_clean %>%
  dplyr::select(
    snp,
    var,
    scaffold,
    position,
    annotated_scaffold,
    gene_id,
    gene_name,
    gene_biotype,
    description,
    gene_start = start,
    gene_end = end,
    strand
  ) %>%
  dplyr::arrange(var, scaffold, position)

lfmm_gene_list <- lfmm_snp_gene_table %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::distinct(gene_id, gene_name, gene_biotype, description) %>%
  dplyr::arrange(gene_name)

write.csv(
  lfmm_snp_gene_table,
  "results/11_lfmm/LFMM_outlier_GENES.csv",
  row.names = FALSE
)



cat("LFMM outlier SNPs:", nrow(lfmm_outliers_mapped), "\n")
cat("SNP-gene overlap rows:", nrow(lfmm_snp_gene_table), "\n")
cat("Unique SNPs inside genes:", length(unique(lfmm_snp_gene_table$snp)), "\n")
cat("Unique genes containing LFMM SNPs:", nrow(lfmm_gene_list), "\n")




































