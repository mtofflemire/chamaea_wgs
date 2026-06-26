#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

output_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/06_Genome-wide_stats"
dataset <- "Chamaea_auto_filteredQC"

fst_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/06_Genome-wide_stats/Chamaea_auto_filteredQC_North_vs_South_10kb.windowed.weir.fst"
north_pi_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/06_Genome-wide_stats/Chamaea_auto_filteredQC_North_10kb.windowed.pi"
south_pi_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/06_Genome-wide_stats/Chamaea_auto_filteredQC_South_10kb.windowed.pi"
north_tajima_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/06_Genome-wide_stats/Chamaea_auto_filteredQC_North_10kb.Tajima.D"
south_tajima_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/06_Genome-wide_stats/Chamaea_auto_filteredQC_South_10kb.Tajima.D"

make_positions <- function(x) {
    x$CHROM <- factor(x$CHROM, levels = unique(x$CHROM))

    chrom_offsets <- x %>%
        group_by(CHROM) %>%
        summarise(chr_offset = max(BIN_START, na.rm = TRUE), .groups = "drop") %>%
        mutate(chr_offset = cumsum(lag(chr_offset, default = 0)))

    x %>%
        left_join(chrom_offsets, by = "CHROM") %>%
        mutate(CUMULATIVE_POS = BIN_START + chr_offset)
}

plot_window_stat <- function(x, y_col, y_label, out_file) {
    threshold_95 <- quantile(x[[y_col]], 0.95, na.rm = TRUE)
    threshold_99 <- quantile(x[[y_col]], 0.99, na.rm = TRUE)

    p <- ggplot(x, aes(x = CUMULATIVE_POS, y = .data[[y_col]], color = CHROM)) +
        geom_point(alpha = 0.5, size = 0.4) +
        geom_hline(yintercept = threshold_95, linetype = "dashed", color = "red") +
        geom_hline(yintercept = threshold_99, linetype = "dashed", color = "black") +
        scale_color_manual(values = rep(c("black", "gray60"), length.out = length(unique(x$CHROM)))) +
        labs(x = "Scaffold", y = y_label) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank()
        )

    ggsave(out_file, p, width = 9, height = 4, dpi = 600)
}

weighted_mean <- function(x, value_col, weight_col) {
    keep <- is.finite(x[[value_col]]) & is.finite(x[[weight_col]]) & x[[weight_col]] > 0
    weighted.mean(x[[value_col]][keep], x[[weight_col]][keep], na.rm = TRUE)
}

fst <- read.table(fst_file, header = TRUE)
north_pi <- read.table(north_pi_file, header = TRUE)
south_pi <- read.table(south_pi_file, header = TRUE)
north_tajima <- read.table(north_tajima_file, header = TRUE)
south_tajima <- read.table(south_tajima_file, header = TRUE)

fst <- make_positions(fst)
north_pi <- make_positions(north_pi)
south_pi <- make_positions(south_pi)
north_tajima <- make_positions(north_tajima)
south_tajima <- make_positions(south_tajima)

north_pi$N_SITES <- north_pi$N_VARIANTS + north_pi$N_MONOMORPHIC
south_pi$N_SITES <- south_pi$N_VARIANTS + south_pi$N_MONOMORPHIC

summary_stats <- data.frame(
    dataset = dataset,
    statistic = c(
        "North_vs_South_Fst",
        "North_pi",
        "South_pi",
        "North_TajimasD",
        "South_TajimasD"
    ),
    weighted_mean = c(
        weighted_mean(fst, "WEIGHTED_FST", "N_VARIANTS"),
        weighted_mean(north_pi, "PI", "N_SITES"),
        weighted_mean(south_pi, "PI", "N_SITES"),
        weighted_mean(north_tajima, "TajimaD", "N_SNPS"),
        weighted_mean(south_tajima, "TajimaD", "N_SNPS")
    ),
    unweighted_mean = c(
        mean(fst$WEIGHTED_FST, na.rm = TRUE),
        mean(north_pi$PI, na.rm = TRUE),
        mean(south_pi$PI, na.rm = TRUE),
        mean(north_tajima$TajimaD, na.rm = TRUE),
        mean(south_tajima$TajimaD, na.rm = TRUE)
    ),
    median = c(
        median(fst$WEIGHTED_FST, na.rm = TRUE),
        median(north_pi$PI, na.rm = TRUE),
        median(south_pi$PI, na.rm = TRUE),
        median(north_tajima$TajimaD, na.rm = TRUE),
        median(south_tajima$TajimaD, na.rm = TRUE)
    )
)

write.csv(
    summary_stats,
    paste0(output_dir, "/", dataset, "_10kb_windowed_stats_summary.csv"),
    row.names = FALSE
)

sink(paste0(output_dir, "/", dataset, "_10kb_windowed_stats_summary.txt"))
print(summary_stats)
sink()

plot_window_stat(
    fst,
    "WEIGHTED_FST",
    "Fst (North vs South)",
    paste0(output_dir, "/", dataset, "_North_vs_South_10kb_Fst_manhattan.png")
)

plot_window_stat(
    north_pi,
    "PI",
    expression(pi~"(North)"),
    paste0(output_dir, "/", dataset, "_North_10kb_pi_manhattan.png")
)

plot_window_stat(
    south_pi,
    "PI",
    expression(pi~"(South)"),
    paste0(output_dir, "/", dataset, "_South_10kb_pi_manhattan.png")
)

plot_window_stat(
    north_tajima,
    "TajimaD",
    "Tajima's D (North)",
    paste0(output_dir, "/", dataset, "_North_10kb_TajimaD_manhattan.png")
)

plot_window_stat(
    south_tajima,
    "TajimaD",
    "Tajima's D (South)",
    paste0(output_dir, "/", dataset, "_South_10kb_TajimaD_manhattan.png")
)
