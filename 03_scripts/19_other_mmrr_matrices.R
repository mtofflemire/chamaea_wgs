# Load necessary library
library(geosphere)  # Provides great-circle distance calculations

# Define file paths
coord_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/chamaea.coords.txt'
geodist_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/geodist.txt'

# Read coordinates (NO HEADERS, JUST LAT/LONG)
coords <- read.table(coord_file, header = FALSE, sep = "\t")

# Compute pairwise great-circle distance matrix (Haversine formula)
geoMat <- distm(coords[,2:3], fun = distHaversine)  # Computes distances in meters

# Save matrix without headers or row names
write.table(geoMat, geodist_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Saved clean great-circle geographic distance matrix:", geodist_file, "\n")














setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea")
# Load libraries
library(raster)

# Your sample order
sample_ids <- readLines("scripts/samples.txt")

# Load metadata and match to sample order
metadata <- read.csv(
  "data/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- metadata[match(sample_ids, metadata$sampleID), ]

# Coordinates in the same order as your genetic matrix
coords_df <- data.frame(
  SampleID = metadata$sampleID,
  Longitude = as.numeric(metadata$long),
  Latitude = as.numeric(metadata$lat)
)

# Load bioclim rasters
bioclim_files <- list.files(
  path = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared/layers/climate/bio_2-5m_bil",
  pattern = "\\.bil$",
  full.names = TRUE
)

bio_nums <- as.numeric(gsub("bio|\\.bil", "", basename(bioclim_files)))
bioclim_files <- bioclim_files[order(bio_nums)]
bio_nums <- bio_nums[order(bio_nums)]

climate <- raster::stack(bioclim_files)
names(climate) <- paste0("bio", bio_nums)

# Extract environmental values at each sample coordinate
pred <- as.data.frame(
  raster::extract(climate, coords_df[, c("Longitude", "Latitude")])
)

colnames(pred) <- names(climate)
rownames(pred) <- coords_df$SampleID

# Split temperature and precipitation variables
temp_pred <- pred[, paste0("bio", 1:11), drop = FALSE]
prec_pred <- pred[, paste0("bio", 12:19), drop = FALSE]

# PCA
pca_temp <- prcomp(temp_pred, center = TRUE, scale. = TRUE)
pca_prec <- prcomp(prec_pred, center = TRUE, scale. = TRUE)

# Keep first 3 PCs from each
temp_pc <- as.data.frame(pca_temp$x[, 1:3])
prec_pc <- as.data.frame(pca_prec$x[, 1:3])

colnames(temp_pc) <- paste0("Temp_PC", 1:3)
colnames(prec_pc) <- paste0("Prec_PC", 1:3)

# This is your environmental PC table
env_pc_fixed <- cbind(temp_pc, prec_pc)
rownames(env_pc_fixed) <- coords_df$SampleID

# THIS is the environmental distance matrix for MMRR
ecoMat <- as.matrix(dist(env_pc_fixed, method = "euclidean"))

# Optional: save it
write.table(
  ecoMat,
  "scripts/ecodist_envPCs.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)






















# Load PLINK PCA eigenvec file
pca <- read.table(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/02_pca/autoAll/Chamaea_autoALL_filteredQC_maf0.05_hwe0.01_prune1kb.pca.eigenvec",
  header = FALSE,
  stringsAsFactors = FALSE
)

# PLINK eigenvec columns are usually:
# V1 = family ID, V2 = sample ID, V3 = PC1, V4 = PC2, etc.
colnames(pca)[1:5] <- c("FID", "SampleID", "PC1", "PC2", "PC3")

# Keep sample IDs and PC1
pc1_scores <- pca[, c("SampleID", "PC1")]

# Optional: if you have a sample order file, match to it
sample_ids <- readLines(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/samples.txt"
)

pc1_scores <- pc1_scores[match(sample_ids, pc1_scores$SampleID), ]

# Check that all samples matched
if (any(is.na(pc1_scores$SampleID))) {
  stop("Some samples in samples.txt were not found in the PCA eigenvec file.")
}

# Make PCA1 distance matrix
pc1Mat <- as.matrix(dist(pc1_scores$PC1, method = "euclidean"))

# Optional names
rownames(pc1Mat) <- pc1_scores$SampleID
colnames(pc1Mat) <- pc1_scores$SampleID

# Optional: save clean matrix for MMRR
write.table(
  pc1Mat,
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/pc1dist.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


















