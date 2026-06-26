required_packages <- c(
  "tidyverse", "gtools", "patchwork", "sf", "cowplot",
  "tigris", "rnaturalearth", "ggforce", "igraph"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg)
  }
}

library(tidyverse)
library(gtools)
library(patchwork)
library(sf)
library(cowplot)
library(tigris)
library(rnaturalearth)
library(ggforce)
library(igraph)

options(tigris_use_cache = TRUE)

args_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(args_file)) dirname(args_file) else getwd()

project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
base_dir <- file.path(project_dir, "4_analyses", "02_admixture")
prefix <- "Chamaea_auto_filteredQC_maf0.05_hwe0.01_prune1kb"

famfile <- file.path(base_dir, paste0(prefix, ".fam"))
metadata_file <- file.path(project_dir, "1_Meta", "CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv")
range_shp <- file.path(project_dir, "2_data", "range_polygon", "SppDataRequest.shp")

k_values <- c(2, 3)

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
  "ancestral3" = "green"
)

fam <- read.table(famfile, header = FALSE)
metadata <- read.csv(metadata_file)

sample_metadata <-
  tibble(sampleID = fam$V2) %>%
  left_join(
    metadata %>%
      select(sampleID, Ecoregion, lat, long),
    by = "sampleID"
  ) %>%
  filter(!is.na(Ecoregion), !is.na(lat), !is.na(long)) %>%
  mutate(
    Ecoregion = factor(Ecoregion, levels = ecoregion_order),
    EcoAbbrev = factor(
      ecoregion_labels[match(Ecoregion, ecoregion_order)],
      levels = ecoregion_labels
    )
  )

# Establish K = 2 sample order and use it to align K = 3 clusters.
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

align_to_k2 <- function(K) {
  qfile <- file.path(base_dir, paste0(prefix, ".", K, ".Q"))
  original_names <- paste0("raw", seq_len(K))

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

make_plot_data <- function(K) {
  cluster_names <- paste0("ancestral", seq_len(K))

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

bar_data <-
  map_dfr(k_values, make_plot_data) %>%
  mutate(
    K_label = factor(K_label, levels = c("K = 2", "K = 3")),
    cluster = factor(as.character(cluster), levels = rev(paste0("ancestral", 1:3)))
  )

p_bars <-
  ggplot(bar_data, aes(sampleID, ancestry, fill = cluster)) +
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
    panel.spacing.y = unit(0.55, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, size = 12, face = "plain"),
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 10)
  )

make_map_data <- function(K) {
  cluster_names <- paste0("ancestral", seq_len(K))

  align_to_k2(K) %>%
    left_join(sample_metadata, by = "sampleID") %>%
    select(sampleID, Ecoregion, EcoAbbrev, long, lat, all_of(cluster_names)) %>%
    mutate(K_label = paste0("K = ", K))
}

map_crs <- 3310
collapse_nearby_samples <- TRUE
cluster_distance_m <- 10000

q2_for_map <-
  make_map_data(2) %>%
  select(sampleID, Ecoregion, EcoAbbrev, long, lat, ancestral1, ancestral2)

q2_sf <-
  q2_for_map %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(map_crs)

cluster_ids <- st_is_within_distance(q2_sf, q2_sf, dist = cluster_distance_m)
cluster_graph <- igraph::graph_from_adj_list(cluster_ids, mode = "all")
q2_sf$locality_cluster <- igraph::components(cluster_graph)$membership

pie_coords <- st_coordinates(q2_sf)

label_offset <- tibble(
  EcoAbbrevShort = ecoregion_abbrev,
  dx = c(
    26000, -26000, -26000, 26000, 26000,
    26000, 30000, 30000, 26000, 28000,
    -30000, -30000, -26000, 30000, 26000
  ),
  dy = c(
    -2000, -6000, 10000, 18000, 16000,
    -12000, 12000, 12000, -14000, -10000,
    26000, -10000, -12000, 10000, -16000
  )
)

main_label_data <-
  q2_sf %>%
  st_drop_geometry() %>%
  mutate(
    x = pie_coords[, 1],
    y = pie_coords[, 2],
    EcoAbbrevShort = ecoregion_abbrev[match(Ecoregion, ecoregion_order)],
    EcoNumber = match(Ecoregion, ecoregion_order)
  ) %>%
  group_by(Ecoregion, EcoAbbrevShort, EcoNumber) %>%
  summarise(
    x = mean(x),
    y = mean(y),
    .groups = "drop"
  ) %>%
  left_join(label_offset, by = "EcoAbbrevShort") %>%
  mutate(
    label = as.character(EcoNumber),
    x_label = x + dx,
    y_label = y + dy
  )

locality_label_data <-
  q2_sf %>%
  st_drop_geometry() %>%
  mutate(
    x = pie_coords[, 1],
    y = pie_coords[, 2],
    EcoAbbrevShort = ecoregion_abbrev[match(Ecoregion, ecoregion_order)]
  ) %>%
  group_by(locality_cluster) %>%
  summarise(
    x = mean(x),
    y = mean(y),
    EcoAbbrevShort = names(sort(table(EcoAbbrevShort), decreasing = TRUE))[1],
    .groups = "drop"
  )

if (collapse_nearby_samples) {
  pie_data <-
    q2_sf %>%
    st_drop_geometry() %>%
    mutate(x = pie_coords[, 1], y = pie_coords[, 2]) %>%
    group_by(locality_cluster) %>%
    summarise(
      x = mean(x),
      y = mean(y),
      n = n(),
      ancestral1 = mean(ancestral1),
      ancestral2 = mean(ancestral2),
      .groups = "drop"
    ) %>%
    mutate(
      pie_id = paste0("locality_", locality_cluster),
      r = 13500 + sqrt(n) * 4700
    ) %>%
    pivot_longer(
      cols = starts_with("ancestral"),
      names_to = "cluster",
      values_to = "ancestry"
    ) %>%
    group_by(pie_id) %>%
    arrange(cluster, .by_group = TRUE) %>%
    mutate(
      end = cumsum(ancestry) * 2 * pi,
      start = lag(end, default = 0),
      r0 = 0
    ) %>%
    ungroup()
} else {
  pie_data <-
    q2_sf %>%
    st_drop_geometry() %>%
    mutate(
      x = pie_coords[, 1],
      y = pie_coords[, 2],
      pie_id = sampleID,
      n = 1,
      r = 17000
    ) %>%
    select(pie_id, x, y, n, r, ancestral1, ancestral2) %>%
    pivot_longer(
      cols = starts_with("ancestral"),
      names_to = "cluster",
      values_to = "ancestry"
    ) %>%
    group_by(pie_id) %>%
    arrange(cluster, .by_group = TRUE) %>%
    mutate(
      end = cumsum(ancestry) * 2 * pi,
      start = lag(end, default = 0),
      r0 = 0
    ) %>%
    ungroup()
}

map_buffer_m <- 130000
map_extent <-
  q2_sf %>%
  st_union() %>%
  st_buffer(dist = map_buffer_m) %>%
  st_as_sf()

map_bbox <- st_bbox(map_extent)
map_bbox["ymin"] <- map_bbox["ymin"] - 90000
map_bbox["ymax"] <- map_bbox["ymax"] + 80000
map_bbox["xmax"] <- map_bbox["xmax"] + 170000

states_sf <-
  tigris::states(cb = TRUE, year = 2022, progress_bar = FALSE) %>%
  filter(!STUSPS %in% c("AK", "HI", "PR")) %>%
  st_transform(map_crs)

countries_sf <-
  rnaturalearth::ne_countries(
    country = c("United States of America", "Canada", "Mexico"),
    scale = "medium",
    returnclass = "sf"
  ) %>%
  st_transform(map_crs)

counties_sf <-
  tigris::counties(state = c("CA", "OR", "NV"), cb = TRUE, year = 2022, progress_bar = FALSE) %>%
  st_transform(map_crs)

range_sf <- NULL
if (file.exists(range_shp)) {
  range_sf <-
    st_read(range_shp, quiet = TRUE) %>%
    filter(SCI_NAME == "Chamaea fasciata") %>%
    st_transform(map_crs) %>%
    st_make_valid() %>%
    st_crop(map_bbox)
}

bay_box_ll <- st_as_sfc(
  st_bbox(
    c(xmin = -123.45, ymin = 36.65, xmax = -120.85, ymax = 39.35),
    crs = st_crs(4326)
  )
) %>%
  st_transform(map_crs)
bay_bbox <- st_bbox(bay_box_ll)
bay_cx <- mean(c(bay_bbox["xmin"], bay_bbox["xmax"]))
bay_cy <- mean(c(bay_bbox["ymin"], bay_bbox["ymax"]))
bay_half_side <- max(
  bay_bbox["xmax"] - bay_bbox["xmin"],
  bay_bbox["ymax"] - bay_bbox["ymin"]
) / 2
bay_bbox["xmin"] <- bay_cx - bay_half_side
bay_bbox["xmax"] <- bay_cx + bay_half_side
bay_bbox["ymin"] <- bay_cy - bay_half_side
bay_bbox["ymax"] <- bay_cy + bay_half_side

inset_label_data <-
  locality_label_data %>%
  filter(
    x >= bay_bbox["xmin"], x <= bay_bbox["xmax"],
    y >= bay_bbox["ymin"], y <= bay_bbox["ymax"]
  ) %>%
  mutate(
    EcoAbbrevShort = if_else(
      EcoAbbrevShort == "NCCR" & y == min(y[EcoAbbrevShort == "NCCR"], na.rm = TRUE),
      "NCICR",
      EcoAbbrevShort
    )
  ) %>%
  mutate(
    label_dx = case_when(
      EcoAbbrevShort == "NCC" ~ -52000,
      EcoAbbrevShort == "CCC" ~ 30000,
      EcoAbbrevShort == "NCICR" ~ 26000,
      EcoAbbrevShort == "NCCR" ~ 26000,
      EcoAbbrevShort == "GV" ~ 30000,
      TRUE ~ 23000
    ),
    label_dy = case_when(
      EcoAbbrevShort == "NCC" ~ -10000,
      EcoAbbrevShort == "CCC" ~ -23000,
      EcoAbbrevShort == "NCICR" ~ -26000,
      EcoAbbrevShort == "NCCR" ~ -26000,
      EcoAbbrevShort == "GV" ~ 23000,
      TRUE ~ 18000
    )
  )

p_map_projected <-
  ggplot() +
  geom_sf(
    data = countries_sf,
    fill = "white",
    color = "black",
    linewidth = 0.50
  )

if (!is.null(range_sf) && nrow(range_sf) > 0) {
  p_map_projected <-
    p_map_projected +
    geom_sf(
      data = range_sf,
      fill = "grey70",
      color = NA,
      alpha = 0.45
    )
}

p_map_projected <-
  p_map_projected +
  geom_rect(
    aes(
      xmin = bay_bbox["xmin"], xmax = bay_bbox["xmax"],
      ymin = bay_bbox["ymin"], ymax = bay_bbox["ymax"]
    ),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
      linewidth = 0.78
  ) +
  geom_sf(
    data = counties_sf,
    fill = NA,
    color = "gray70",
    linewidth = 0.12
  ) +
  geom_sf(
    data = states_sf,
    fill = NA,
    color = "black",
    linewidth = 0.46
  ) +
  geom_sf(
    data = countries_sf,
    fill = NA,
    color = "black",
    linewidth = 0.70
  ) +
  geom_arc_bar(
    data = pie_data,
    aes(
      x0 = x,
      y0 = y,
      r0 = r0,
      r = r,
      start = start,
      end = end,
      fill = cluster
    ),
    color = "black",
    linewidth = 0.34,
    alpha = 0.98
  ) +
  geom_text(
    data = main_label_data,
    aes(x = x_label, y = y_label, label = label),
    inherit.aes = FALSE,
    size = 3.0,
    color = "black"
  ) +
  scale_fill_manual(
    values = cluster_palette[c("ancestral1", "ancestral2")],
    labels = c("Cluster 1", "Cluster 2"),
    name = "K = 2"
  ) +
  coord_sf(
    xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
    ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
    crs = st_crs(map_crs),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_map_inset <-
  ggplot() +
  geom_sf(
    data = countries_sf,
    fill = "white",
    color = "black",
    linewidth = 0.70
  )

if (!is.null(range_sf) && nrow(range_sf) > 0) {
  p_map_inset <-
    p_map_inset +
    geom_sf(
      data = range_sf,
      fill = "grey70",
      color = NA,
      alpha = 0.45
    )
}

pie_data_inset <-
  pie_data %>%
  mutate(r = r * 0.80)

p_map_inset <-
  p_map_inset +
  geom_sf(
    data = counties_sf,
    fill = NA,
    color = "gray70",
    linewidth = 0.14
  ) +
  geom_sf(
    data = states_sf,
    fill = NA,
    color = "black",
    linewidth = 0.45
  ) +
  geom_sf(
    data = countries_sf,
    fill = NA,
    color = "black",
    linewidth = 0.70
  ) +
  geom_arc_bar(
    data = pie_data_inset,
    aes(
      x0 = x,
      y0 = y,
      r0 = r0,
      r = r,
      start = start,
      end = end,
      fill = cluster
    ),
    color = "black",
    linewidth = 0.34,
    alpha = 0.98
  ) +
  geom_text(
    data = inset_label_data,
    aes(
      x = x + label_dx,
      y = y + label_dy,
      label = EcoAbbrevShort
    ),
    inherit.aes = FALSE,
    size = 2.7,
    color = "black"
  ) +
  scale_fill_manual(
    values = cluster_palette[c("ancestral1", "ancestral2")],
    labels = c("Cluster 1", "Cluster 2"),
    name = "K = 2"
  ) +
  coord_sf(
    xlim = c(bay_bbox["xmin"], bay_bbox["xmax"]),
    ylim = c(bay_bbox["ymin"], bay_bbox["ymax"]),
    crs = st_crs(map_crs),
    expand = FALSE
  ) +
  annotate(
    "text",
    x = bay_bbox["xmin"] + 25000,
    y = bay_bbox["ymin"] + 25000,
    label = "San Francisco\nBay Area",
    hjust = 0,
    vjust = 0,
    size = 3.0
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1.00)
  )

p_map_projected_with_inset <-
  ggdraw() +
  draw_plot(p_map_projected, x = 0.34, y = 0.00, width = 0.62, height = 1.00) +
  draw_plot(p_map_inset, x = 0.025, y = 0.60, width = 0.32, height = 0.32) +
  theme(plot.background = element_rect(fill = "white", color = NA))

bar_pdf <- file.path(base_dir, paste0(prefix, "_ADMIXTURE_K2-K3_by_ecoregion_right_labels.pdf"))
bar_png <- file.path(base_dir, paste0(prefix, "_ADMIXTURE_K2-K3_by_ecoregion_right_labels.png"))
map_pdf <- file.path(base_dir, paste0(prefix, "_ADMIXTURE_K2_map_projected_collapsed_with_inset.pdf"))
map_png <- file.path(base_dir, paste0(prefix, "_ADMIXTURE_K2_map_projected_collapsed_with_inset.png"))
inset_pdf <- file.path(base_dir, paste0(prefix, "_ADMIXTURE_K2_map_BayArea_inset_only.pdf"))
inset_png <- file.path(base_dir, paste0(prefix, "_ADMIXTURE_K2_map_BayArea_inset_only.png"))

ggsave(bar_pdf, p_bars, width = 10.5, height = 3.6, useDingbats = FALSE)
ggsave(bar_png, p_bars, width = 10.5, height = 3.6, dpi = 600)

ggsave(map_pdf, p_map_projected_with_inset, width = 7.6, height = 7.0, useDingbats = FALSE)
ggsave(map_png, p_map_projected_with_inset, width = 7.6, height = 7.0, dpi = 600)
ggsave(inset_pdf, p_map_inset, width = 4.2, height = 4.2, useDingbats = FALSE)
ggsave(inset_png, p_map_inset, width = 4.2, height = 4.2, dpi = 600)

message("Saved: ", bar_pdf)
message("Saved: ", bar_png)
message("Saved: ", map_pdf)
message("Saved: ", map_png)
message("Saved: ", inset_pdf)
message("Saved: ", inset_png)
