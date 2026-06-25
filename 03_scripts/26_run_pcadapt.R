args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript 03_run_pcadapt.R ROOT")
root <- normalizePath(args[1], mustWork = TRUE)

suppressPackageStartupMessages({
  library(pcadapt)
  library(qvalue)
  library(data.table)
  library(ggplot2)
})

prefix <- file.path(root, "results", "pcadapt", "Chamaea_confirmed_prune1kb")
bed_file <- paste0(prefix, ".bed")
coords_file <- file.path(root, "results", "pcadapt", "variant_coordinates.tsv")
out_dir <- file.path(root, "results", "pcadapt")
log_file <- file.path(root, "logs", "03_pcadapt.log")

sink(log_file, split = TRUE)
cat("PCAdapt version:", as.character(packageVersion("pcadapt")), "\n")
cat("Input:", bed_file, "\n")
cat("K: 1\n")

genotype_data <- read.pcadapt(bed_file, type = "bed")
fit_file <- file.path(out_dir, "pcadapt_K1_fit.rds")
if (file.exists(fit_file)) {
  cat("Loading existing PCAdapt fit:", fit_file, "\n")
  fit <- readRDS(fit_file)
} else {
  fit <- pcadapt(genotype_data, K = 1)
  saveRDS(fit, fit_file, compress = FALSE)
}

coords <- fread(coords_file)
if (nrow(coords) != length(fit$pvalues)) {
  stop("Coordinate rows do not match PCAdapt p-value count.")
}

qvals <- qvalue(fit$pvalues)$qvalues
bonf <- p.adjust(fit$pvalues, method = "bonferroni")
bh <- p.adjust(fit$pvalues, method = "BH")

stats <- data.table(
  SNP_INDEX = seq_along(fit$pvalues),
  coords,
  P_VALUE = fit$pvalues,
  Q_VALUE = qvals,
  BONFERRONI_P = bonf,
  BH_P = bh,
  MAF = fit$maf,
  PASS = seq_along(fit$pvalues) %in% fit$pass
)
fwrite(stats, file.path(out_dir, "pcadapt_K1_all_statistics.tsv"), sep = "\t")

bonf_out <- stats[!is.na(BONFERRONI_P) & BONFERRONI_P < 0.05]
q_out <- stats[!is.na(Q_VALUE) & Q_VALUE < 0.01]
bh_out <- stats[!is.na(BH_P) & BH_P < 0.01]
fwrite(bonf_out, file.path(out_dir, "pcadapt_K1_bonferroni_0.05_outliers.tsv"), sep = "\t")
fwrite(q_out, file.path(out_dir, "pcadapt_K1_qvalue_0.01_outliers.tsv"), sep = "\t")
fwrite(bh_out, file.path(out_dir, "pcadapt_K1_BH_0.01_outliers.tsv"), sep = "\t")

summary <- data.table(
  metric = c("input_variants", "tested_variants", "bonferroni_0.05_outliers",
             "qvalue_0.01_outliers", "BH_0.01_outliers", "genomic_inflation_factor"),
  value = c(nrow(stats), sum(!is.na(stats$P_VALUE)), nrow(bonf_out),
            nrow(q_out), nrow(bh_out), fit$gif)
)
fwrite(summary, file.path(out_dir, "pcadapt_K1_summary.tsv"), sep = "\t")

plot_dt <- copy(stats)
plot_dt[, LOGP := -log10(P_VALUE)]
plot_dt[, scaffold_order := match(CHROM, unique(CHROM))]
lengths <- plot_dt[, .(max_pos = max(POS)), by = .(CHROM, scaffold_order)][order(scaffold_order)]
lengths[, offset := shift(cumsum(max_pos), fill = 0)]
plot_dt <- merge(plot_dt, lengths[, .(CHROM, offset)], by = "CHROM", sort = FALSE)
plot_dt[, CUM_POS := POS + offset]
plot_dt[, color_group := scaffold_order %% 2]

p <- ggplot(plot_dt, aes(CUM_POS, LOGP, color = factor(color_group))) +
  geom_point(size = 0.25, alpha = 0.65) +
  geom_hline(yintercept = -log10(0.05 / sum(!is.na(stats$P_VALUE))),
             color = "#D62728", linewidth = 0.35) +
  scale_color_manual(values = c("#00A83B", "#8DDE00")) +
  labs(x = NULL, y = expression(-log[10](italic(p))),
       title = "PCAdapt K = 1 (Bonferroni threshold)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(file.path(out_dir, "pcadapt_K1_manhattan.png"), p,
       width = 11, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, "pcadapt_K1_manhattan.pdf"), p,
       width = 11, height = 4.5)

png(file.path(out_dir, "pcadapt_K1_diagnostics.png"),
    width = 2400, height = 1800, res = 250)
par(mfrow = c(2, 2))
plot(fit, option = "scores")
plot(fit, option = "qqplot")
hist(fit$pvalues, breaks = 50, xlab = "p-values", main = "P-value distribution")
plot(fit, option = "stat.distribution")
dev.off()

cat("Bonferroni outliers:", nrow(bonf_out), "\n")
cat("Q-value outliers:", nrow(q_out), "\n")
sink()
