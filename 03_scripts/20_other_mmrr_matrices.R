library(geosphere)
library(raster)

sample_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/3_scripts/samples.txt"
metadata_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/1_Meta/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv"
pca_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/01_pca/Chamaea_auto_filteredQC_maf0.05_hwe0.01_prune1kb.pca.eigenvec"
climate_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared/layers/climate/wc2"
out_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/09_mmrr"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- readLines(sample_file)

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- metadata[match(samples, metadata$sampleID), ]

if (any(is.na(metadata$sampleID))) {
  stop("Some samples were not found in the metadata file.")
}

coords <- data.frame(
  sampleID = metadata$sampleID,
  long = as.numeric(metadata$long),
  lat = as.numeric(metadata$lat)
)

write.table(coords, paste0(out_dir, "/mmrr_sample_coordinates.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

geoMat <- distm(coords[, c("long", "lat")], fun = distHaversine)
rownames(geoMat) <- coords$sampleID
colnames(geoMat) <- coords$sampleID

climate_files <- paste0(climate_dir, "/wc2.1_2.5m_bio_", 1:19, ".tif")
if (!all(file.exists(climate_files))) {
  stop("Some WorldClim raster files are missing.")
}

climate <- stack(climate_files)
names(climate) <- paste0("bio", 1:19)

env <- as.data.frame(raster::extract(climate, coords[, c("long", "lat")]))
rownames(env) <- coords$sampleID

if (any(is.na(env))) {
  stop("Missing climate values were extracted for at least one sample.")
}

temp_pca <- prcomp(env[, paste0("bio", 1:11)], center = TRUE, scale. = TRUE)
prec_pca <- prcomp(env[, paste0("bio", 12:19)], center = TRUE, scale. = TRUE)

env_pcs <- data.frame(
  sampleID = coords$sampleID,
  Temp_PC1 = temp_pca$x[, 1],
  Temp_PC2 = temp_pca$x[, 2],
  Temp_PC3 = temp_pca$x[, 3],
  Prec_PC1 = prec_pca$x[, 1],
  Prec_PC2 = prec_pca$x[, 2],
  Prec_PC3 = prec_pca$x[, 3]
)

write.table(env_pcs, paste0(out_dir, "/mmrr_climate_pc_scores.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

ecoMat <- as.matrix(dist(env_pcs[, 2:7], method = "euclidean"))
rownames(ecoMat) <- coords$sampleID
colnames(ecoMat) <- coords$sampleID

pca <- read.table(pca_file, header = TRUE, stringsAsFactors = FALSE)
pca <- pca[match(samples, pca$IID), ]

if (any(is.na(pca$IID))) {
  stop("Some samples were not found in the PCA eigenvec file.")
}

pc1Mat <- as.matrix(dist(pca$PC1, method = "euclidean"))
rownames(pc1Mat) <- pca$IID
colnames(pc1Mat) <- pca$IID

genMat <- as.matrix(read.table(paste0(out_dir, "/gendist.all.txt"), header = FALSE))
rownames(genMat) <- samples
colnames(genMat) <- samples

write_matrix <- function(x, file) {
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

writeLines(samples, paste0(out_dir, "/samples.all.txt"))
write_matrix(genMat, paste0(out_dir, "/gendist.all.txt"))
write_matrix(geoMat, paste0(out_dir, "/geodist.all.txt"))
write_matrix(ecoMat, paste0(out_dir, "/ecodist.All.txt"))
write_matrix(pc1Mat, paste0(out_dir, "/pc1dist.all.txt"))

sink(paste0(out_dir, "/mmrr_matrix_summary.txt"))
cat("All samples:", length(samples), "\n")
cat("Temperature PCA variance explained:\n")
print(summary(temp_pca)$importance[, 1:3])
cat("\nPrecipitation PCA variance explained:\n")
print(summary(prec_pca)$importance[, 1:3])
sink()
