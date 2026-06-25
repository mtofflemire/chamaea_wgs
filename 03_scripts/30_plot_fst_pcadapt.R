args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 09_make_four_panel_selection_figure.R ROOT")
}

root <- normalizePath(args[1], mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggrastr)
  library(patchwork)
  library(jsonlite)
})

pc_file <- file.path(
  root, "results", "pcadapt", "pcadapt_K1_all_statistics.tsv"
)
ns_file <- file.path(
  root, "results", "fst",
  "North_vs_South_50kb_step10kb.windowed.weir.fst"
)
or_nccr_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_results/Selective_sweeps/OC_vs_NCCR_fst.windowed.weir.fst"
or_snf_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_results/Selective_sweeps/OC_vs_SNF_fst.windowed.weir.fst"
gene_file <- file.path(
  root, "inputs", "reference_genes_mapped_to_vcf_scaffolds.tsv"
)
sequence_report <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/2_data/reference/ncbi_dataset/data/GCF_029207755.1/sequence_report.jsonl"

out_dir <- file.path(root, "results", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c("#E41A1C", "#FC8D62", "#8DD321", "#2C7FB8")
label_genes <- c(
  "TPCN2", "NADSYN1", "DHCR7", "MORC2", "RAB3GAP2",
  "SMARCA4", "PPFIBP2", "LUZP2", "DPP7", "FGD5",
  "CDKN1B", "LOC136290045"
)

pc <- fread(pc_file)
pc <- pc[!is.na(P_VALUE) & is.finite(P_VALUE)]

sequence_lines <- readLines(sequence_report)
sequence_lengths <- rbindlist(lapply(sequence_lines, function(line) {
  item <- fromJSON(line)
  data.table(CHROM = item$genbankAccession, scaffold_length = item$length)
}), fill = TRUE)
sequence_lengths <- unique(sequence_lengths[!is.na(CHROM)])

scaffold_order <- unique(pc$CHROM)
scaffolds <- data.table(CHROM = scaffold_order, scaffold_index = seq_along(scaffold_order))
scaffolds <- merge(scaffolds, sequence_lengths, by = "CHROM", all.x = TRUE, sort = FALSE)
setorder(scaffolds, scaffold_index)
scaffolds[is.na(scaffold_length), scaffold_length := 0]
scaffolds[, offset := shift(cumsum(scaffold_length), fill = 0)]
scaffolds[, color_group := factor((scaffold_index - 1L) %% 4L)]

add_coordinates <- function(data, position_column) {
  result <- merge(
    data,
    scaffolds[, .(CHROM, scaffold_index, offset, color_group)],
    by = "CHROM",
    all.x = TRUE,
    sort = FALSE
  )
  result[, CUM_POS := get(position_column) + offset]
  result
}

pc <- add_coordinates(pc, "POS")
pc[, LOGP := -log10(P_VALUE)]
pc[, QVALUE_OUTLIER := !is.na(Q_VALUE) & Q_VALUE < 0.01]
q_boundary <- max(pc[QVALUE_OUTLIER == TRUE, P_VALUE])
pc_line <- -log10(q_boundary)

read_fst <- function(path) {
  data <- fread(path)
  data <- data[is.finite(WEIGHTED_FST)]
  data[, MIDPOINT := (as.numeric(BIN_START) + as.numeric(BIN_END)) / 2]
  add_coordinates(data, "MIDPOINT")
}

ns <- read_fst(ns_file)
or_nccr <- read_fst(or_nccr_file)
or_snf <- read_fst(or_snf_file)

fst_threshold <- function(data) {
  mean(data$WEIGHTED_FST) + 5 * sd(data$WEIGHTED_FST)
}

ns_threshold <- fst_threshold(ns)
or_nccr_threshold <- fst_threshold(or_nccr)
or_snf_threshold <- fst_threshold(or_snf)

genes <- fread(gene_file)
genes <- genes[GENE_SYMBOL %in% label_genes]

pc_gene_hits <- merge(
  pc[, .(CHROM, POS, CUM_POS, LOGP, Q_VALUE)],
  genes[, .(CHROM, START, END, GENE_SYMBOL)],
  by = "CHROM",
  allow.cartesian = TRUE
)[POS >= START & POS <= END & Q_VALUE < 0.01]
gene_fst_labels <- function(fst_data) {
  hits <- merge(
    fst_data[, .(
      CHROM, BIN_START, BIN_END, MIDPOINT, CUM_POS, WEIGHTED_FST
    )],
    genes[, .(CHROM, START, END, GENE_SYMBOL)],
    by = "CHROM",
    allow.cartesian = TRUE
  )[
    BIN_START <= END & BIN_END >= START
  ]
  hits[, .SD[which.max(WEIGHTED_FST)], by = GENE_SYMBOL]
}

ns_gene_hits <- merge(
  ns[, .(
    CHROM, BIN_START, BIN_END, MIDPOINT, CUM_POS, WEIGHTED_FST
  )],
  genes[, .(CHROM, START, END, GENE_SYMBOL)],
  by = "CHROM",
  allow.cartesian = TRUE
)[BIN_START <= END & BIN_END >= START]

left_peak_groups <- data.table(
  PEAK_GROUP = c(
    "RAB3GAP2", "CDKN1B_CLUSTER", "TPCN2_CLUSTER", "FGD5",
    "MORC2", "DPP7", "SMARCA4"
  ),
  LABEL = c(
    "RAB3GAP2",
    "CDKN1B\nLOC136290045",
    "DHCR7\nPPFIBP2\nNADSYN1",
    "FGD5",
    "MORC2",
    "DPP7",
    "SMARCA4"
  ),
  GENES = c(
    "RAB3GAP2",
    "CDKN1B;LOC136290045",
    "DHCR7;PPFIBP2;NADSYN1",
    "FGD5",
    "MORC2",
    "DPP7",
    "SMARCA4"
  ),
  X_SHIFT_BP = c(
    0, 0, 0, 0, -6000000, 7000000, 0
  )
)

make_pc_peak_labels <- function() {
  rbindlist(lapply(seq_len(nrow(left_peak_groups)), function(i) {
    gene_set <- strsplit(left_peak_groups$GENES[i], ";", fixed = TRUE)[[1]]
    hits <- if (left_peak_groups$PEAK_GROUP[i] == "TPCN2_CLUSTER") {
      pc[
        CHROM == "JARCOQ010000007.1" &
          POS >= 8270000 & POS <= 8720000 &
          Q_VALUE < 0.01
      ]
    } else if (left_peak_groups$PEAK_GROUP[i] == "CDKN1B_CLUSTER") {
      pc[
        CHROM == "JARCOQ010000005.1" &
          POS >= 74350000 & POS <= 74450000 &
          Q_VALUE < 0.01
      ]
    } else {
      pc_gene_hits[GENE_SYMBOL %in% gene_set]
    }
    if (!nrow(hits)) return(NULL)
    peak <- hits[which.max(LOGP)]
    data.table(
      PEAK_GROUP = left_peak_groups$PEAK_GROUP[i],
      LABEL = left_peak_groups$LABEL[i],
      CUM_POS = peak$CUM_POS,
      LABEL_X = peak$CUM_POS + left_peak_groups$X_SHIFT_BP[i],
      PEAK_Y = peak$LOGP,
      LABEL_Y = peak$LOGP + 0.55
    )
  }))
}

make_ns_peak_labels <- function() {
  rbindlist(lapply(seq_len(nrow(left_peak_groups)), function(i) {
    is_cluster <- left_peak_groups$PEAK_GROUP[i] == "TPCN2_CLUSTER"
    gene_set <- if (is_cluster) {
      c("TPCN2", "NADSYN1", "DHCR7")
    } else {
      strsplit(left_peak_groups$GENES[i], ";", fixed = TRUE)[[1]]
    }
    hits <- ns_gene_hits[GENE_SYMBOL %in% gene_set]
    if (!nrow(hits)) return(NULL)
    peak <- hits[which.max(WEIGHTED_FST)]
    data.table(
      PEAK_GROUP = left_peak_groups$PEAK_GROUP[i],
      LABEL = if (is_cluster) {
        "TPCN2\nNADSYN1\nDHCR7"
      } else {
        left_peak_groups$LABEL[i]
      },
      CUM_POS = peak$CUM_POS,
      LABEL_X = peak$CUM_POS + left_peak_groups$X_SHIFT_BP[i],
      PEAK_Y = peak$WEIGHTED_FST,
      LABEL_Y = min(0.84, peak$WEIGHTED_FST + 0.025)
    )
  }))
}

pc_labels <- make_pc_peak_labels()
ns_labels <- make_ns_peak_labels()

fst_outlier_gene_hits <- function(fst_data, threshold) {
  outliers <- fst_data[WEIGHTED_FST >= threshold]
  hits <- merge(
    outliers[, .(
      CHROM, BIN_START, BIN_END, MIDPOINT, CUM_POS, WEIGHTED_FST
    )],
    fread(gene_file)[, .(
      CHROM, START, END, GENE_ID, GENE_SYMBOL, DESCRIPTION
    )],
    by = "CHROM",
    allow.cartesian = TRUE
  )[BIN_START <= END & BIN_END >= START]
  hits
}

nccr_gene_hits <- fst_outlier_gene_hits(or_nccr, or_nccr_threshold)
snf_gene_hits <- fst_outlier_gene_hits(or_snf, or_snf_threshold)

nccr_gene_summary <- nccr_gene_hits[, .(
  MAX_NCCR_FST = max(WEIGHTED_FST),
  N_NCCR_WINDOWS = uniqueN(paste(CHROM, BIN_START, BIN_END))
), by = .(GENE_ID, GENE_SYMBOL, DESCRIPTION)]
snf_gene_summary <- snf_gene_hits[, .(
  MAX_SNF_FST = max(WEIGHTED_FST),
  N_SNF_WINDOWS = uniqueN(paste(CHROM, BIN_START, BIN_END))
), by = .(GENE_ID, GENE_SYMBOL, DESCRIPTION)]

shared_or_genes <- merge(
  nccr_gene_summary,
  snf_gene_summary,
  by = c("GENE_ID", "GENE_SYMBOL", "DESCRIPTION")
)
shared_or_genes[, MIN_SHARED_PEAK := pmin(MAX_NCCR_FST, MAX_SNF_FST)]
shared_or_genes[, MEAN_SHARED_PEAK := (MAX_NCCR_FST + MAX_SNF_FST) / 2]
setorder(shared_or_genes, -MIN_SHARED_PEAK, -MEAN_SHARED_PEAK)

# Stack named genes by shared genomic divergence peak.
or_peak_groups <- data.table(
  PEAK_GROUP = c("G3BP2_USO1", "TTN"),
  LABEL = c("G3BP2\nUSO1", "TTN"),
  GENES = c("G3BP2;USO1", "TTN"),
  X_SHIFT_BP = c(0, 18000000)
)

make_or_peak_labels <- function(gene_hits) {
  labels <- rbindlist(lapply(seq_len(nrow(or_peak_groups)), function(i) {
    gene_set <- strsplit(or_peak_groups$GENES[i], ";", fixed = TRUE)[[1]]
    hits <- gene_hits[GENE_SYMBOL %in% gene_set]
    if (!nrow(hits)) return(NULL)
    peak <- hits[which.max(WEIGHTED_FST)]
    data.table(
      PEAK_GROUP = or_peak_groups$PEAK_GROUP[i],
      LABEL = or_peak_groups$LABEL[i],
      CUM_POS = peak$CUM_POS,
      PEAK_Y = peak$WEIGHTED_FST,
      LABEL_X = peak$CUM_POS + or_peak_groups$X_SHIFT_BP[i],
      LABEL_Y = min(0.975, peak$WEIGHTED_FST + 0.045)
    )
  }))
  labels
}

or_nccr_labels <- make_or_peak_labels(nccr_gene_hits)
or_snf_labels <- make_or_peak_labels(snf_gene_hits)
shared_or_label_genes <- unique(unlist(strsplit(
  or_peak_groups$GENES, ";", fixed = TRUE
)))

base_theme <- theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(linewidth = 0.45),
    plot.margin = margin(7, 8, 7, 8)
  )

point_layer <- function(data, y_column, point_size = 0.28) {
  geom_point_rast(
    data = data,
    aes(x = CUM_POS, y = .data[[y_column]], color = color_group),
    size = point_size,
    alpha = 0.72,
    raster.dpi = 600
  )
}

panel_a <- ggplot() +
  point_layer(pc, "LOGP", 0.45) +
  geom_hline(
    yintercept = pc_line, color = "#E41A1C",
    linewidth = 0.35
  ) +
  geom_segment(
    data = pc_labels[!grepl("\n", LABEL, fixed = TRUE)],
    aes(
      x = LABEL_X, xend = CUM_POS,
      y = LABEL_Y - 0.12, yend = PEAK_Y + 0.08
    ),
    linewidth = 0.3,
    color = "gray25"
  ) +
  geom_text(
    data = pc_labels,
    aes(LABEL_X, LABEL_Y, label = LABEL),
    size = 2.65,
    fontface = "italic",
    color = "black",
    lineheight = 0.82,
    vjust = 0
  ) +
  scale_color_manual(values = palette) +
  labs(y = expression(-log[10](italic(p)) ~ "(PCAdapt)"), tag = "A") +
  base_theme

panel_b <- ggplot() +
  point_layer(ns, "WEIGHTED_FST", 0.50) +
  geom_hline(
    yintercept = ns_threshold, color = "#E41A1C",
    linewidth = 0.35
  ) +
  geom_segment(
    data = ns_labels[!grepl("\n", LABEL, fixed = TRUE)],
    aes(
      x = LABEL_X, xend = CUM_POS,
      y = LABEL_Y - 0.006, yend = PEAK_Y + 0.004
    ),
    linewidth = 0.3,
    color = "gray25"
  ) +
  geom_text(
    data = ns_labels,
    aes(LABEL_X, LABEL_Y, label = LABEL),
    size = 2.65,
    fontface = "italic",
    color = "black",
    lineheight = 0.82,
    vjust = 0
  ) +
  scale_color_manual(values = palette) +
  coord_cartesian(ylim = c(0, max(ns$WEIGHTED_FST) * 1.08)) +
  labs(y = "Fst (North vs South)", tag = "B") +
  base_theme

panel_c <- ggplot() +
  point_layer(or_nccr, "WEIGHTED_FST", 0.50) +
  geom_hline(
    yintercept = or_nccr_threshold, color = "#E41A1C",
    linewidth = 0.35
  ) +
  geom_segment(
    data = or_nccr_labels[!grepl("\n", LABEL, fixed = TRUE)],
    aes(
      x = LABEL_X, xend = CUM_POS,
      y = LABEL_Y - 0.006, yend = PEAK_Y + 0.004
    ),
    linewidth = 0.3,
    color = "gray25"
  ) +
  geom_text(
    data = or_nccr_labels,
    aes(LABEL_X, LABEL_Y, label = LABEL),
    size = 2.65,
    fontface = "italic",
    color = "black",
    lineheight = 0.82,
    vjust = 0
  ) +
  scale_color_manual(values = palette) +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(y = "Fst (OR vs NCCR)", tag = "C") +
  base_theme

panel_d <- ggplot() +
  point_layer(or_snf, "WEIGHTED_FST", 0.50) +
  geom_hline(
    yintercept = or_snf_threshold, color = "#E41A1C",
    linewidth = 0.35
  ) +
  geom_segment(
    data = or_snf_labels[!grepl("\n", LABEL, fixed = TRUE)],
    aes(
      x = LABEL_X, xend = CUM_POS,
      y = LABEL_Y - 0.006, yend = PEAK_Y + 0.004
    ),
    linewidth = 0.3,
    color = "gray25"
  ) +
  geom_text(
    data = or_snf_labels,
    aes(LABEL_X, LABEL_Y, label = LABEL),
    size = 2.65,
    fontface = "italic",
    color = "black",
    lineheight = 0.82,
    vjust = 0
  ) +
  scale_color_manual(values = palette) +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(y = "Fst (OR vs SNF)", tag = "D") +
  base_theme

combined <- (panel_a | panel_c) / (panel_b | panel_d) &
  theme(
    plot.tag = element_text(face = "bold", size = 19),
    plot.tag.position = c(0.01, 0.98)
  )

png_file <- file.path(
  out_dir,
  "updated_PCadapt_Fst_four_panel_selection_figure.png"
)
pdf_file <- file.path(
  out_dir,
  "updated_PCadapt_Fst_four_panel_selection_figure.pdf"
)
hybrid_pdf_file <- file.path(
  out_dir,
  "updated_PCadapt_Fst_four_panel_selection_figure_hybrid_600dpi.pdf"
)

ggsave(
  png_file, combined,
  width = 16, height = 7.4, dpi = 400,
  bg = "white"
)
ggsave(
  pdf_file, combined,
  width = 16, height = 7.4,
  bg = "white"
)
ggsave(
  hybrid_pdf_file, combined,
  device = cairo_pdf,
  width = 16, height = 7.4,
  bg = "white"
)

thresholds <- data.table(
  panel = c("A", "B", "C", "D"),
  analysis = c(
    "PCAdapt K=1 q-value < 0.01",
    "North vs South Fst",
    "Oregon Coast vs NCCR Fst",
    "Oregon Coast vs SNF Fst"
  ),
  threshold = c(
    pc_line, ns_threshold, or_nccr_threshold, or_snf_threshold
  ),
  threshold_scale = c(
    "-log10(raw p-value boundary corresponding to q < 0.01)",
    rep("mean weighted Fst + 5 SD", 3)
  )
)
fwrite(
  thresholds,
  file.path(out_dir, "updated_four_panel_figure_thresholds.tsv"),
  sep = "\t"
)

fwrite(
  data.table(GENE_SYMBOL = label_genes),
  file.path(out_dir, "updated_four_panel_labeled_genes.tsv"),
  sep = "\t"
)

fwrite(
  shared_or_genes,
  file.path(out_dir, "genes_shared_between_OR_NCCR_and_OR_SNF_outlier_windows.tsv"),
  sep = "\t"
)
fwrite(
  data.table(GENE_SYMBOL = shared_or_label_genes),
  file.path(out_dir, "genes_labeled_in_both_OR_Fst_panels.tsv"),
  sep = "\t"
)

cat("Saved:", png_file, "\n")
cat("Saved:", pdf_file, "\n")
cat("Saved:", hybrid_pdf_file, "\n")
