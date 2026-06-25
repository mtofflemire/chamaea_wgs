# load libraries
library(tidyverse)
library(gtools)
library(patchwork)

# file paths
base_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_results/2_admixture/run1_seed43"
prefix <- "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb"

famfile <- file.path(base_dir, paste0(prefix, ".fam"))
metadata_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/1_Meta/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv"

k_values <- 2:10
stacked_k_values <- 4:10



ecoregion_order <- c(
  "Oregon Coast", "Northern California Coast", "Klamath Mountains",
  "Northern California Coast Ranges", "Northern California Interior Coast Ranges",
  "Southern Cascades", "Great Valley", "Modoc Plateau",
  "Sierra Nevada", "Sierra Nevada Foothills", "Central Valley Coast Ranges",
  "Central California Coast", "Southern California Coast",
  "Southern California Mountains and Valleys", "Colorado Desert"
)



ecoregion_abbrev <- c(
  "OR", "NCC", "KM", "NCCR", "NCICR", "SC", "GV", "MP",
  "SN", "SNF", "CVCR", "CCC", "SCC", "SCMV", "CD"
)



ecoregion_labels <- paste0(ecoregion_abbrev, " [", seq_along(ecoregion_abbrev), "]")



cluster_palette <- c(
  "ancestral1" = "blue",
  "ancestral2" = "red",
  "ancestral3" = "green",
  "ancestral4" = "purple",
  "ancestral5" = "orange",
  "ancestral6" = "cyan3",
  "ancestral7" = "magenta",
  "ancestral8" = "gold",
  "ancestral9" = "gray30",
  "ancestral10" = "brown"
)

# load shared data
fam <- read.table(famfile, header = FALSE)
metadata <- read.csv(metadata_file)

sample_metadata <-
  tibble(sampleID = fam$V2) %>%
  left_join(metadata %>% select(sampleID, Ecoregion), by = "sampleID") %>%
  filter(!is.na(Ecoregion)) %>%
  mutate(
    Ecoregion = factor(Ecoregion, levels = ecoregion_order),
    EcoAbbrev = factor(
      ecoregion_labels[match(Ecoregion, ecoregion_order)],
      levels = ecoregion_labels
    )
  )

# establish K = 2 sample order and cluster identity
q2_raw <-
  read.table(file.path(base_dir, paste0(prefix, ".2.Q")), header = FALSE) %>%
  as_tibble() %>%
  setNames(c("ancestral1", "ancestral2")) %>%
  mutate(sampleID = fam$V2) %>%
  left_join(sample_metadata, by = "sampleID") %>%
  filter(!is.na(EcoAbbrev)) %>%
  arrange(EcoAbbrev, desc(ancestral1))

sample_order <- q2_raw$sampleID

q2_matrix <-
  q2_raw %>%
  select(sampleID, ancestral1, ancestral2) %>%
  arrange(match(sampleID, sample_order))

# align K clusters to K = 2 clusters
align_to_k2 <- function(K) {
  
  qfile <- file.path(base_dir, paste0(prefix, ".", K, ".Q"))
  original_names <- paste0("raw", 1:K)
  
  qk_raw <-
    read.table(qfile, header = FALSE) %>%
    as_tibble() %>%
    setNames(original_names) %>%
    mutate(sampleID = fam$V2) %>%
    filter(sampleID %in% sample_order) %>%
    arrange(match(sampleID, sample_order))
  
  if (K == 2) {
    return(
      qk_raw %>%
        transmute(sampleID, ancestral1 = raw1, ancestral2 = raw2)
    )
  }
  
  q2_values <- as.matrix(q2_matrix[, c("ancestral1", "ancestral2")])
  qk_values <- as.matrix(qk_raw[, original_names])
  
  perms <- permutations(n = K, r = 2, v = seq_len(K))
  
  scores <- apply(perms, 1, function(p) {
    cor(q2_values[, 1], qk_values[, p[1]]) +
      cor(q2_values[, 2], qk_values[, p[2]])
  })
  
  best <- perms[which.max(scores), ]
  extra <- setdiff(seq_len(K), best)
  
  ordered_cols <- c(best, extra)
  aligned_names <- paste0("ancestral", seq_len(K))
  
  qk_raw %>%
    select(sampleID, all_of(original_names[ordered_cols])) %>%
    setNames(c("sampleID", aligned_names))
}

# prepare plotting data
make_plot_data <- function(K) {
  
  cluster_names <- paste0("ancestral", 1:K)
  
  align_to_k2(K) %>%
    left_join(sample_metadata, by = "sampleID") %>%
    mutate(
      sampleID = factor(sampleID, levels = sample_order),
      K_label = factor(paste0("K = ", K), levels = paste0("K = ", k_values))
    ) %>%
    arrange(sampleID) %>%
    pivot_longer(
      cols = all_of(cluster_names),
      names_to = "cluster",
      values_to = "ancestry"
    ) %>%
    mutate(cluster = factor(cluster, levels = rev(cluster_names)))
}

# individual plot function
plot_admixture_individual <- function(K) {
  
  cluster_names <- paste0("ancestral", 1:K)
  plot_data <- make_plot_data(K)
  
  ggplot(plot_data, aes(sampleID, ancestry, fill = cluster)) +
    geom_col(color = "gray30", linewidth = 0.08, width = 1) +
    facet_grid(~ EcoAbbrev, switch = "x", scales = "free_x", space = "free_x") +
    labs(title = paste0("K = ", K), y = "Ancestry", x = NULL) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.7)) +
    scale_fill_manual(values = cluster_palette[cluster_names], guide = "none") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 10, face = "bold"),
      plot.title.position = "plot",
      panel.spacing.x = unit(0.08, "lines"),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7)
    )
}

# Plot K = 2 and K = 3 together, with K labels on the right.
plot_data_k2_k3 <-
  map_dfr(c(2, 3), make_plot_data) %>%
  mutate(
    K_label = factor(K_label, levels = c("K = 2", "K = 3")),
    cluster = factor(as.character(cluster),
                     levels = rev(paste0("ancestral", 1:3)))
  )

p <-
  ggplot(plot_data_k2_k3, aes(sampleID, ancestry, fill = cluster)) +
  geom_col(color = "gray30", linewidth = 0.08, width = 1) +
  facet_grid(
    K_label ~ EcoAbbrev,
    switch = "x",
    scales = "free_x",
    space = "free_x"
  ) +
  labs(y = "Ancestry", x = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.7)) +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  theme_minimal() +
  theme(
    panel.spacing.x = unit(0.08, "lines"),
    panel.spacing.y = unit(0.12, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5, size = 7
    ),
    strip.text.y.right = element_text(
      angle = 0, hjust = 0.5, size = 12, face = "plain"
    ),
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 10)
  )

pdf_file <- file.path(
  base_dir,
  paste0(prefix, "_ADMIXTURE_K2-K3_by_ecoregion_right_labels.pdf")
)
png_file <- file.path(
  base_dir,
  paste0(prefix, "_ADMIXTURE_K2-K3_by_ecoregion_right_labels.png")
)

ggsave(pdf_file, p, width = 10.5, height = 4.0, useDingbats = FALSE)
ggsave(png_file, p, width = 10.5, height = 4.0, dpi = 600)

message("Saved: ", pdf_file)
message("Saved: ", png_file)
quit(save = "no", status = 0)


















library(ggplot2)

log_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/out/2_ADMIXTURE/autoAll"

log_files <- list.files(
  log_dir,
  pattern = "^log[0-9]+\\.out$",
  full.names = TRUE
)

parse_cv <- function(file) {
  x <- readLines(file)
  cv_line <- grep("^CV error", x, value = TRUE)
  
  data.frame(
    K = as.numeric(sub(".*K=([0-9]+).*", "\\1", cv_line)),
    CV_error = as.numeric(sub(".*: *([0-9.]+).*", "\\1", cv_line))
  )
}

cv_results <- do.call(rbind, lapply(log_files, parse_cv))
cv_results <- cv_results[order(cv_results$K), ]

p <- ggplot(cv_results, aes(x = K, y = CV_error)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  labs(
    x = "K",
    y = "CV error"
  ) +
  theme_bw()

ggsave(
  filename = file.path(log_dir, "admixture_cv_error_plot.pdf"),
  plot = p,
  width = 7,
  height = 5
)

ggsave(
  filename = file.path(log_dir, "admixture_cv_error_plot.png"),
  plot = p,
  width = 7,
  height = 4,
  dpi = 300
)








# 4x2 patchwork layout for K4-K10
admixture_plots <- map(stacked_k_values, plot_admixture_individual)

combined_plot <-
  wrap_plots(
    admixture_plots,
    ncol = 2,
    byrow = FALSE
  )

combined_plot

ggsave(
  file.path(base_dir, paste0(prefix, "_ADMIXTURE_K4-K10_4x2_K2_aligned.pdf")),
  combined_plot,
  width = 20,
  height = 12,
  useDingbats = FALSE
)

ggsave(
  file.path(base_dir, paste0(prefix, "_ADMIXTURE_K4-K10_4x2_K2_aligned.png")),
  combined_plot,
  width = 8,
  height = 9,
  dpi = 600
)
