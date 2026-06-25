# load libraries
required_packages <- c("ggplot2", "data.table", "scales", "scatterplot3d",
                       "maps", "akima", "dplyr", "sf")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}



# file paths
command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
script_dir <- if (length(script_arg)) {
  script_path <- sub("^--file=", "", script_arg[1])
  script_path <- gsub("~\\+~", " ", script_path)
  dirname(normalizePath(script_path))
} else {
  normalizePath(getwd())
}

project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
base_dir <- file.path(project_dir, "4_analyses", "01_pca")

pca_eigenvec_file <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb.pca.eigenvec")
pca_eigenval_file <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb.pca.eigenval")
metadata_file <- file.path(
  project_dir, "1_Meta",
  "CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv"
)
range_file <- file.path(
  project_dir, "2_data", "SppDataRequest", "SppDataRequest.shp"
)

output_3d_pdf <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb_PCA_3D_plot.pdf")
output_3d_png <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb_PCA_3D_plot.png")
output_interp_pdf <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb_PC1_interpolated_map.pdf")
output_interp_png <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb_PC1_interpolated_map.png")

# load data
pca_data <- fread(pca_eigenvec_file, header = TRUE)
eigenvalues <- fread(pca_eigenval_file, header = FALSE)
meta_data <- fread(metadata_file, header = TRUE)

merged_data <- merge(pca_data, meta_data, by.x = "IID", by.y = "sampleID", all.x = TRUE)
merged_data <- merged_data[complete.cases(merged_data[, c("PC1", "PC2", "PC3", "long", "lat", "Ecoregion")]), ]
merged_data$Ecoregion <- factor(merged_data$Ecoregion)

# PCA labels and colors
percent_variance <- (eigenvalues$V1 / sum(eigenvalues$V1)) * 100

x_label <- paste0("PC1 (", round(percent_variance[1], 2), "% variance)")
y_label <- paste0("PC2 (", round(percent_variance[2], 2), "% variance)")
z_label <- paste0("PC3 (", round(percent_variance[3], 2), "% variance)")

ecoregion_colours <- c(
  "black", "#E69F00", "#56B4E9", "green", "#F0E442",
  "#0072B2", "purple4", "hotpink", "#999999", "blue",
  "orangered", "#4DAF4A", "#984EA3", "red3", "#FFFF33"
)

ecoregion_colours <- alpha(ecoregion_colours, 0.7)

# 3D PCA plot
plot_pca_3d <- function(data, legend_inset = c(0.055, 0.045)) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  ecoregion_levels <- levels(data$Ecoregion)
  colour_lookup <- setNames(ecoregion_colours[seq_along(ecoregion_levels)], ecoregion_levels)
  point_colours <- colour_lookup[as.character(data$Ecoregion)]
  
  par(mar = c(4.6, 4.6, 1.2, 3.2), xpd = FALSE)
  
  scatterplot3d(
    x = data$PC1,
    y = data$PC2,
    z = data$PC3,
    color = point_colours,
    pch = 19,
    xlab = x_label,
    ylab = "",
    zlab = z_label,
    angle = 58,
    scale.y = 0.9,
    y.axis.offset = 0.75,
    type = "h",
    box = TRUE,
    cex.symbols = 0.8,
    cex.axis = 0.8,
    cex.lab = 0.9
  )

  mtext(y_label, side = 4, line = -0.8, cex = 0.9)
  
  legend(
    "topright",
    inset = legend_inset,
    legend = ecoregion_levels,
    col = colour_lookup,
    pch = 19,
    cex = 0.48,
    pt.cex = 0.7,
    bty = "o",
    bg = "white",
    box.col = "black",
    box.lwd = 1,
    x.intersp = 0.65,
    y.intersp = 1.12
  )
}

pdf(output_3d_pdf, width = 9.3, height = 4.95)
plot_pca_3d(merged_data, legend_inset = c(0.10, 0.005))
dev.off()

png(output_3d_png, width = 3441, height = 1832, res = 300)
plot_pca_3d(merged_data)
dev.off()

# Load the Wrentit range polygon. The source shapefile contains several
# species, so select Chamaea fasciata explicitly.
species_ranges <- st_read(range_file, quiet = TRUE)
wrentit_range <- species_ranges[
  species_ranges$SCI_NAME == "Chamaea fasciata",
]

if (!nrow(wrentit_range)) {
  stop("Chamaea fasciata was not found in: ", range_file)
}

wrentit_range <- st_make_valid(wrentit_range)
wrentit_range <- st_transform(wrentit_range, 4326)

# Interpolate PC1 only where supported by the sampling locations. Cells
# outside the sampling convex hull remain blank, avoiding unsupported
# extrapolation across unsampled portions of the species range.
interp_result <- with(
  merged_data,
  interp(
    x = long, y = lat, z = PC1,
    duplicate = "mean",
    linear = TRUE,
    extrap = FALSE,
    nx = 300,
    ny = 300
  )
)

interp_data <- expand.grid(
  Longitude = interp_result$x,
  Latitude = interp_result$y
)
interp_data$PC1 <- as.vector(interp_result$z)
interp_data <- interp_data[!is.na(interp_data$PC1), ]

interp_points <- st_as_sf(
  interp_data,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)
inside_wrentit_range <- lengths(st_intersects(
  interp_points, wrentit_range
)) > 0
interp_data <- interp_data[inside_wrentit_range, ]

# interpolated PC1 map
plot_pc1_map <-
  ggplot() +
  borders("world", regions = c("USA", "Mexico"),
          fill = "gray95", colour = "gray20", linewidth = 0.3) +
  borders("state", colour = "gray20", linewidth = 0.15) +
  geom_raster(data = interp_data,
              aes(x = Longitude, y = Latitude, fill = PC1),
              alpha = 0.8) +
  scale_fill_gradientn(colors = c("blue", "green", "red"), name = "PC1") +
  coord_sf(xlim = c(-125, -113), ylim = c(29, 46.5), expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

ggsave(output_interp_pdf, plot = plot_pc1_map, width = 5, height = 7)
ggsave(output_interp_png, plot = plot_pc1_map, width = 5, height = 7, dpi = 600)

print(plot_pc1_map)
























#now plot for supscpies colors and add to supplemental
# load libraries
required_packages <- c("ggplot2", "data.table", "scales", "scatterplot3d",
                       "maps", "akima", "dplyr", "sf")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# file paths
base_dir <- file.path(project_dir, "4_analyses", "01_pca")

pca_eigenvec_file <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb.pca.eigenvec")
pca_eigenval_file <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb.pca.eigenval")
metadata_file <- file.path(
  project_dir, "1_Meta",
  "CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv"
)

output_3d_pdf <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb_PCA_3D_plot_subspecies.pdf")
output_3d_png <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb_PCA_3D_plot_subspecies.png")
output_interp_pdf <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb_PC1_interpolated_map.pdf")
output_interp_png <- file.path(base_dir, "Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb_PC1_interpolated_map.png")

# load data
pca_data <- fread(pca_eigenvec_file, header = TRUE)
eigenvalues <- fread(pca_eigenval_file, header = FALSE)
meta_data <- fread(metadata_file, header = TRUE)

merged_data <- merge(pca_data, meta_data, by.x = "IID", by.y = "sampleID", all.x = TRUE)

merged_data <- merged_data[
  complete.cases(
    merged_data[, c("PC1", "PC2", "PC3", "long", "lat", "Ecoregion", "subspecies")]
  ),
]

merged_data$Ecoregion <- factor(merged_data$Ecoregion)
merged_data$subspecies <- factor(merged_data$subspecies)

# PCA labels and colors
percent_variance <- (eigenvalues$V1 / sum(eigenvalues$V1)) * 100

x_label <- paste0("PC1 (", round(percent_variance[1], 2), "% variance)")
y_label <- paste0("PC2 (", round(percent_variance[2], 2), "% variance)")
z_label <- paste0("PC3 (", round(percent_variance[3], 2), "% variance)")

subspecies_colours <- c(
  "#000000",  # black
  "#D7191C",  # strong red
  "#2C7BB6",  # strong blue
  "#1A9641",  # strong green
  "#762A83",  # dark purple
  "#Fdae61"   # orange
)

subspecies_colours <- alpha(subspecies_colours, 1)


subspecies_colours <- alpha(subspecies_colours, 0.75)

# 3D PCA plot colored by subspecies
plot_pca_3d <- function(data, legend_inset = c(0.055, 0.045)) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  subspecies_levels <- levels(data$subspecies)
  colour_lookup <- setNames(
    subspecies_colours[seq_along(subspecies_levels)],
    subspecies_levels
  )
  
  point_colours <- colour_lookup[as.character(data$subspecies)]
  
  par(mar = c(4.6, 4.6, 1.2, 3.2), xpd = FALSE)
  
  scatterplot3d(
    x = data$PC1,
    y = data$PC2,
    z = data$PC3,
    color = point_colours,
    pch = 19,
    xlab = x_label,
    ylab = "",
    zlab = z_label,
    angle = 58,
    scale.y = 0.9,
    y.axis.offset = 0.75,
    type = "h",
    box = TRUE,
    cex.symbols = 0.8,
    cex.axis = 0.8,
    cex.lab = 0.9
  )

  mtext(y_label, side = 4, line = -0.8, cex = 0.9)
  
  legend(
    "topright",
    inset = legend_inset,
    legend = subspecies_levels,
    col = colour_lookup,
    pch = 19,
    cex = 0.52,
    pt.cex = 0.7,
    bty = "o",
    bg = "white",
    box.col = "black",
    box.lwd = 1,
    x.intersp = 0.65,
    y.intersp = 1.12
  )
}

pdf(output_3d_pdf, width = 10.3, height = 5.45)
plot_pca_3d(merged_data, legend_inset = c(0.10, 0.005))
dev.off()

png(output_3d_png, width = 3430, height = 1815, res = 300)
plot_pca_3d(merged_data)
dev.off()
