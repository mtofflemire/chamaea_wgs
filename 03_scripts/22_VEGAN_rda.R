# load libraries
library(raster)
library(vegan)
library(algatr)
library(ggplot2)
library(ggrepel)
library(MASS)
library(dplyr)
library(qvalue)



# set working directory
setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea")


# load samples
sample_ids <- readLines("scripts/samples.txt")


# load coords
metadata <- read.csv("data/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv", stringsAsFactors = FALSE, check.names = FALSE)
metadata <- metadata[match(sample_ids, metadata$sampleID), ]
coords_df <- data.frame(SampleID = metadata$sampleID, Longitude = as.numeric(metadata$long), Latitude = as.numeric(metadata$lat))



# load and process environmental variables
bioclim_files <- list.files(path = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared/layers/climate/bio_2-5m_bil", pattern = "\\.bil$", full.names = TRUE)
bio_nums <- as.numeric(gsub("bio|\\.bil", "", basename(bioclim_files)))
bioclim_files <- bioclim_files[order(bio_nums)]
bio_nums <- bio_nums[order(bio_nums)]
climate <- raster::stack(bioclim_files)
names(climate) <- paste0("bio", bio_nums)



# extract climate values at sample coordinates
pred <- as.data.frame(raster::extract(climate, coords_df[, c("Longitude", "Latitude")]))
colnames(pred) <- names(climate)



# split temperature and precipitation variables
temp_pred <- pred[, paste0("bio", 1:11), drop = FALSE]
prec_pred <- pred[, paste0("bio", 12:19), drop = FALSE]


# PCA
pca_temp <- prcomp(temp_pred, center = TRUE, scale. = TRUE)
pca_prec <- prcomp(prec_pred, center = TRUE, scale. = TRUE)


# retain first 3 PCs for each set
temp_pc <- as.data.frame(pca_temp$x[, 1:3])
prec_pc <- as.data.frame(pca_prec$x[, 1:3])
colnames(temp_pc) <- paste0("Temp_PC", 1:3)
colnames(prec_pc) <- paste0("Prec_PC", 1:3)
env_pc_fixed <- cbind(temp_pc, prec_pc)
rownames(env_pc_fixed) <- coords_df$SampleID
env_pc <- cbind(SampleID = coords_df$SampleID, env_pc_fixed)
env_pc_with_coords <- cbind(coords_df, env_pc_fixed)














#############################################################
### Climate PCA variable contribution plots
### No manual arrows; factoextra-style loading plots
#############################################################
#
#library(factoextra)
#library(ggplot2)
#library(patchwork)
#
#dir.create("results/10_rda/climate_pca_loadings", recursive = TRUE, showWarnings = FALSE)
#
#############################################################
### Temperature PCA
#############################################################
#
#temp_pc12 <- fviz_pca_var(
#  pca_temp,
#  axes = c(1, 2),
#  col.var = "contrib",
#  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#  repel = TRUE
#) +
#  labs(title = "Temperature: PC1 vs PC2") +
#  theme_bw() +
#  theme(
#    plot.title = element_text(hjust = 0.5),
#    panel.grid.minor = element_blank()
#  )
#
#temp_pc13 <- fviz_pca_var(
#  pca_temp,
#  axes = c(1, 3),
#  col.var = "contrib",
#  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#  repel = TRUE
#) +
#  labs(title = "Temperature: PC1 vs PC3") +
#  theme_bw() +
#  theme(
#    plot.title = element_text(hjust = 0.5),
#    panel.grid.minor = element_blank()
#  )
#
#temp_pc23 <- fviz_pca_var(
#  pca_temp,
#  axes = c(2, 3),
#  col.var = "contrib",
#  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#  repel = TRUE
#) +
#  labs(title = "Temperature: PC2 vs PC3") +
#  theme_bw() +
#  theme(
#    plot.title = element_text(hjust = 0.5),
#    panel.grid.minor = element_blank()
#  )
#
#############################################################
### Precipitation PCA
#############################################################
#
#prec_pc12 <- fviz_pca_var(
#  pca_prec,
#  axes = c(1, 2),
#  col.var = "contrib",
#  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#  repel = TRUE
#) +
#  labs(title = "Precipitation: PC1 vs PC2") +
#  theme_bw() +
#  theme(
#    plot.title = element_text(hjust = 0.5),
#    panel.grid.minor = element_blank()
#  )
#
#prec_pc13 <- fviz_pca_var(
#  pca_prec,
#  axes = c(1, 3),
#  col.var = "contrib",
#  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#  repel = TRUE
#) +
#  labs(title = "Precipitation: PC1 vs PC3") +
#  theme_bw() +
#  theme(
#    plot.title = element_text(hjust = 0.5),
#    panel.grid.minor = element_blank()
#  )
#
#prec_pc23 <- fviz_pca_var(
#  pca_prec,
#  axes = c(2, 3),
#  col.var = "contrib",
#  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#  repel = TRUE
#) +
#  labs(title = "Precipitation: PC2 vs PC3") +
#  theme_bw() +
#  theme(
#    plot.title = element_text(hjust = 0.5),
#    panel.grid.minor = element_blank()
#  )
#
#############################################################
### Combine all six plots: 3 rows x 2 columns
#############################################################
#
#climate_pca_loading_plot <- (temp_pc12 | prec_pc12) /
#  (temp_pc13 | prec_pc13) /
#  (temp_pc23 | prec_pc23)
#
#print(climate_pca_loading_plot)
#
#ggsave(
#  "results/10_rda/climate_pca_loadings/climate_PCA_variable_contributions_3x2.pdf",
#  climate_pca_loading_plot,
#  width = 8,
#  height = 10
#)
#
#ggsave(
#  "results/10_rda/climate_pca_loadings/climate_PCA_variable_contributions_3x2.png",
#  climate_pca_loading_plot,
#  width = 8,
#  height = 10,
#  dpi = 600
#)
#







############################################################
## Combined climate PCA RGB composite
## Red = Climate PC1
## Green = Climate PC2
## Blue = Climate PC3
############################################################

library(raster)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)

dir.create("results/10_rda/climate_pca_rgb_maps", recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Run one PCA on all 19 BIOCLIM variables together
############################################################

climate_vars <- paste0("bio", 1:19)

pred_all_climate <- pred[, climate_vars, drop = FALSE]

pca_climate_all <- prcomp(
  pred_all_climate,
  center = TRUE,
  scale. = TRUE
)

print(summary(pca_climate_all))

############################################################
## 2. Project first 3 climate PCs onto raster landscape
############################################################

climate <- climate[[climate_vars]]
names(climate) <- climate_vars

climate_pc_rasters <- raster::predict(
  climate,
  pca_climate_all,
  index = 1:3
)

names(climate_pc_rasters) <- c("Climate_PC1", "Climate_PC2", "Climate_PC3")

############################################################
## 3. Load Wrentit range and mask rasters
############################################################

bird_data <- sf::st_read(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared/shapefiles/SppDataRequest/SppDataRequest.shp",
  quiet = TRUE
)

bird_data_filtered <- bird_data[bird_data$SCI_NAME == "Chamaea fasciata", ]
bird_data_filtered <- sf::st_transform(bird_data_filtered, crs = raster::crs(climate_pc_rasters)@projargs)

bird_range_sp <- as(bird_data_filtered, "Spatial")

climate_pc_crop <- raster::crop(climate_pc_rasters, bird_range_sp)
climate_pc_mask <- raster::mask(climate_pc_crop, bird_range_sp)

############################################################
## 4. Convert PC rasters to RGB values
############################################################

rgb_df <- as.data.frame(climate_pc_mask, xy = TRUE, na.rm = TRUE)

stretch_channel <- function(x, lower_q = 0.02, upper_q = 0.98) {
  qs <- quantile(x, probs = c(lower_q, upper_q), na.rm = TRUE)
  x <- pmin(pmax(x, qs[1]), qs[2])
  scales::rescale(x, to = c(0, 1))
}

rgb_df$R <- stretch_channel(rgb_df$Climate_PC1)
rgb_df$G <- stretch_channel(rgb_df$Climate_PC2)
rgb_df$B <- stretch_channel(rgb_df$Climate_PC3)

rgb_df$RGB <- rgb(rgb_df$R, rgb_df$G, rgb_df$B)

cat("RGB pixels:", nrow(rgb_df), "\n")

############################################################
## 5. Plot RGB composite
############################################################

states <- rnaturalearth::ne_states(
  country = "United States of America",
  returnclass = "sf"
)

climate_rgb_map <- ggplot() +
  geom_raster(
    data = rgb_df,
    aes(x = x, y = y, fill = RGB)
  ) +
  scale_fill_identity() +
  geom_sf(
    data = states,
    fill = NA,
    color = "gray35",
    linewidth = 0.15
  ) +
  geom_sf(
    data = bird_data_filtered,
    fill = NA,
    color = "black",
    linewidth = 0.45
  ) +
  coord_sf(
    xlim = c(-130, -115),
    ylim = c(30, 47),
    expand = FALSE
  ) +
  theme_void() +
  labs(
    title = "Climate PCA RGB Composite"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16)
  )

print(climate_rgb_map)

ggsave(
  "results/10_rda/climate_pca_rgb_maps/combined_climate_PCA_RGB_composite.pdf",
  climate_rgb_map,
  width = 6,
  height = 7
)

ggsave(
  "results/10_rda/climate_pca_rgb_maps/combined_climate_PCA_RGB_composite.png",
  climate_rgb_map,
  width = 6,
  height = 7,
  dpi = 600
)




















#load SNPs
library(vcfR)

# load VCF
vcf <- read.vcfR("results/01_VARIANTS/GATK_QC/AUTOSOMES/Chamaea_autoALL_filteredQC_maf0.05_miss0.25_prune1kb.vcf.gz")

# extract GT strings
gt_raw <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# convert GT strings to allele counts
# 0/0 = 0, 0/1 = 1, 1/1 = 2
gt_raw[gt_raw %in% c(".", "./.", ".|.")] <- NA
gt_raw[gt_raw %in% c("0/0", "0|0")] <- 0
gt_raw[gt_raw %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
gt_raw[gt_raw %in% c("1/1", "1|1")] <- 2

genotypes <- gt_raw
storage.mode(genotypes) <- "numeric"

# VCF gives SNPs x samples; RDA needs samples x SNPs
genotypes <- t(genotypes)

# force same order as scripts/samples.txt
genotypes <- genotypes[sample_ids, , drop = FALSE]

# impute missing values by SNP mode
genotypes <- apply(genotypes, 2, function(x) {replace(x, is.na(x), as.numeric(names(which.max(table(x)))))})
genotypes <- as.matrix(genotypes)
rownames(genotypes) <- sample_ids


















#run RDA
results_full <- rda_run(
  gen = genotypes, env = env_pc_fixed,
  coords = coords_df[, c("Longitude", "Latitude")],
  model = "full", correctGEO = FALSE, correctPC = TRUE, nPC = 2)








#summaryize results
summary(results_full)
saveRDS(results_full,"results/10_rda/rda_full_model_all_loci.rds")
RsquareAdj(results_full)
summary(results_full)$concont
screeplot(results_full)
vif.cca(results_full)
plot(results_full, scaling = 3)
rda_model_test <- anova(results_full, permutations = 999)
rda_axis_test <- anova(results_full, by = "axis", permutations = 999)
rda_term_test <- anova(results_full, by = "term", permutations = 999)
print(rda_model_test)
print(rda_axis_test)
print(rda_term_test)











#now load rda for future use
#setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea")
#results_full <- readRDS("results/10_rda/rda_full_model_all_loci.rds")









#########
#plot rda
#########

# ecoregion colors matching PCA script
ecoregion_colours <- c(
  "black", 
  "#E69F00", 
  "#56B4E9", 
  "green", 
  "#F0E442",
  "#0072B2", 
  "purple4", 
  "hotpink", 
  "#999999", 
  "blue",
  "orangered", 
  "#4DAF4A", 
  "#984EA3", 
  "red3", 
  "#FFFF33"
)



ecoregion_levels <- levels(factor(metadata$Ecoregion))
ecoregion_palette <- setNames(scales::alpha(ecoregion_colours[seq_along(ecoregion_levels)], 0.85),ecoregion_levels)



# sample scores
site_scores <- as.data.frame(scores(results_full, display = "sites", choices = c(1, 2), scaling = 3))
site_scores$SampleID <- rownames(site_scores)
site_scores$Longitude <- coords_df$Longitude
site_scores$Latitude <- coords_df$Latitude
site_scores$Ecoregion <- factor(metadata$Ecoregion, levels = names(ecoregion_palette))



# environmental loading arrows
env_vectors <- as.data.frame(scores(results_full, display = "bp", choices = c(1, 2), scaling = 3))
env_vectors$Variable <- rownames(env_vectors)



# rescale arrows to fit sample points
sample_radius <- min(max(abs(site_scores$RDA1), na.rm = TRUE), max(abs(site_scores$RDA2), na.rm = TRUE)) * 0.75
env_radius <- max(sqrt(env_vectors$RDA1^2 + env_vectors$RDA2^2), na.rm = TRUE)
arrow_multiplier <- sample_radius / env_radius



# plot
rda_ecoregion_arrows_plot <- ggplot() +
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, color = Ecoregion), size = 3, alpha = 0.85) +
  geom_segment(data = env_vectors, aes(x = 0, y = 0, xend = RDA1 * arrow_multiplier, yend = RDA2 * arrow_multiplier), arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.7) +
  geom_text_repel(data = env_vectors, aes(x = RDA1 * arrow_multiplier, y = RDA2 * arrow_multiplier, label = Variable), color = "black", size = 4, max.overlaps = Inf) +
  scale_color_manual(values = ecoregion_palette, drop = FALSE) +
  theme_bw() +
  labs( x = "RDA1", y = "RDA2", color = "Ecoregion") +
  theme(
    legend.position = c(0.6, 1.04),
    legend.justification = c(1, 1),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 2.2, alpha = 1)
    )
  )

print(rda_ecoregion_arrows_plot)




#ggsave("results/10_rda/rda_sampling_points_ecoregion_with_env_arrows.pdf",rda_ecoregion_arrows_plot,width = 5.5,height = 5.5)
#ggsave("results/10_rda/rda_sampling_points_ecoregion_with_env_arrows.png",rda_ecoregion_arrows_plot,width = 5.5,height = 5.5,dpi = 600)









# SNP scores/loadings
snp_scores <- as.data.frame(scores(results_full, display = "species", choices = c(1, 2), scaling = 3))
snp_scores$SNP <- rownames(snp_scores)

# expand SNP cloud visually without moving sample points
snp_zoom <- 7
snp_scores$RDA1_zoom <- snp_scores$RDA1 * snp_zoom
snp_scores$RDA2_zoom <- snp_scores$RDA2 * snp_zoom

# make arrows longer
arrow_multiplier_long <- arrow_multiplier * 2

rda_ecoregion_arrows_plot <- ggplot() +
  geom_point(data = snp_scores,aes(x = RDA1_zoom, y = RDA2_zoom),color = "gray30",alpha = 0.18,size = 0.35) +
  geom_hline(yintercept = 0,linetype = "dashed",color = "gray45",linewidth = 0.35) +
  geom_vline(xintercept = 0,linetype = "dashed",color = "gray45",linewidth = 0.35) +
  geom_point(data = site_scores,aes(x = RDA1, y = RDA2, fill = Ecoregion),shape = 21,color = "gray35",stroke = 0.35,size = 3,alpha = 0.85) +
  geom_segment(data = env_vectors,aes(x = 0,y = 0,xend = RDA1 * arrow_multiplier_long,yend = RDA2 * arrow_multiplier_long),arrow = arrow(length = unit(0.2, "cm")),color = "black",linewidth = 0.7) +
  geom_text_repel(data = env_vectors,aes(x = RDA1 * arrow_multiplier_long,y = RDA2 * arrow_multiplier_long,label = Variable),color = "black",size = 4,max.overlaps = Inf) +
  scale_fill_manual(values = ecoregion_palette, drop = FALSE) +
  labs(x = "RDA1",y = "RDA2",fill = "Ecoregion") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = c(0.4, 0.95),
    legend.justification = c(1, 1),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(ncol = 1,override.aes = list(size = 2.2, alpha = 1)))



print(rda_ecoregion_arrows_plot)



ggsave("results/10_rda/RDA_plot_1.pdf",rda_ecoregion_arrows_plot,width = 8,height = 6)
ggsave("results/10_rda/RDA_plot_1.png",rda_ecoregion_arrows_plot,width = 8,height = 6,dpi = 600)
























##################################
#malahanobis snp outlier detection
##################################
K <- 2
fdr_threshold <- 0.01



# SNP loadings from first K constrained RDA axes
snp_loadings <- as.data.frame(scores(results_full, display = "species", choices = 1:K, scaling = 3))
colnames(snp_loadings) <- paste0("RDA", 1:K)
snp_loadings$SNP <- rownames(snp_loadings)



# z-transform SNP loadings
snp_z <- as.data.frame(scale(snp_loadings[, paste0("RDA", 1:K)]))
colnames(snp_z) <- paste0("RDA", 1:K, "_z")


# robust Mahalanobis distances
set.seed(123)
robust_cov <- MASS::cov.rob(snp_z)
mahal_distances <- mahalanobis(x = snp_z,center = robust_cov$center,cov = robust_cov$cov)



# p-values from chi-square distribution
p_values <- pchisq(mahal_distances,df = K,lower.tail = FALSE)


# q-values / FDR
q_values <- qvalue(p_values)$qvalues


# full results table
mahal_df <- cbind(snp_loadings,snp_z,Mahalanobis_Dist = mahal_distances,p_value = p_values,q_value = q_values)
mahal_df$Chromosome <- sub("_[^_]+$", "", mahal_df$SNP)
mahal_df$Position <- as.numeric(sub("^.*_", "", mahal_df$SNP))



# strict outlier SNPs
outlier_snps <- mahal_df %>%
  filter(q_value < fdr_threshold) %>%
  arrange(q_value, p_value)



# Unique RDA outlier SNPs overall
unique_outlier_snps <- outlier_snps %>%
  dplyr::distinct(SNP, .keep_all = TRUE) %>%
  dplyr::arrange(q_value, p_value)




cat("Number of RDA outlier SNPs at q < 1e-4:", nrow(outlier_snps), "\n")
cat("Total unique RDA outlier SNPs:", nrow(unique_outlier_snps), "\n")



write.csv(mahal_df,"results/10_rda/RDA_results.csv",row.names = FALSE)
write.csv(outlier_snps,"results/10_rda/RDA_outlier_SNPs.csv",row.names = FALSE)
#write.csv(unique_outlier_snps,"results/10_rda/RDA_outlier_SNPs_unique.csv",row.names = FALSE)

























############################################################
## Investigate RDA candidate SNPs
## 1. Strongest RDA loading axis
## 2. Strongest environmental predictor correlation
############################################################

# Use unique SNPs only
cand <- unique_outlier_snps

############################################################
## Count strongest RDA axis/loading for each SNP
############################################################

z_cols <- paste0("RDA", 1:K, "_z")

cand <- cand %>%
  dplyr::mutate(
    strongest_RDA_axis = paste0(
      "RDA",
      max.col(abs(dplyr::select(., dplyr::all_of(z_cols))), ties.method = "first")
    ),
    strongest_RDA_loading = apply(
      abs(dplyr::select(., dplyr::all_of(z_cols))),
      1,
      max
    )
  )

cat("Outlier SNPs by strongest RDA axis:\n")
print(table(cand$strongest_RDA_axis))


############################################################
## Correlate each outlier SNP with environmental PC predictors
############################################################

env_predictors <- as.data.frame(env_pc_fixed)

# Make sure row order matches genotype matrix
env_predictors <- env_predictors[rownames(genotypes), , drop = FALSE]

predictor_names <- colnames(env_predictors)

# Empty matrix for SNP-environment correlations
cor_mat <- matrix(
  NA,
  nrow = nrow(cand),
  ncol = length(predictor_names)
)

colnames(cor_mat) <- predictor_names
rownames(cor_mat) <- cand$SNP

for (i in seq_len(nrow(cand))) {
  snp_id <- cand$SNP[i]
  
  snp_genotype <- genotypes[, snp_id]
  
  for (j in seq_along(predictor_names)) {
    pred_name <- predictor_names[j]
    
    cor_mat[i, j] <- cor(
      snp_genotype,
      env_predictors[[pred_name]],
      use = "complete.obs",
      method = "pearson"
    )
  }
}

cor_df <- as.data.frame(cor_mat)
cor_df$SNP <- rownames(cor_df)

############################################################
## Find strongest environmental predictor for each SNP
############################################################

cand <- cand %>%
  dplyr::left_join(cor_df, by = "SNP")

cand$predictor <- predictor_names[
  max.col(abs(cand[, predictor_names, drop = FALSE]), ties.method = "first")
]

cand$correlation <- apply(
  abs(cand[, predictor_names, drop = FALSE]),
  1,
  max
)

cat("Outlier SNPs by strongest environmental predictor:\n")
print(table(cand$predictor))


############################################################
## Save outputs
############################################################

write.csv(
  cand,
  "results/10_rda/RDA_outlier_SNPs_with_axis_and_predictor.csv",
  row.names = FALSE
)

axis_counts <- as.data.frame(table(cand$strongest_RDA_axis))
colnames(axis_counts) <- c("RDA_axis", "num_outlier_snps")

predictor_counts <- as.data.frame(table(cand$predictor))
colnames(predictor_counts) <- c("predictor", "num_outlier_snps")

write.csv(
  axis_counts,
  "results/10_rda/RDA_outlier_SNP_counts_by_axis.csv",
  row.names = FALSE
)

write.csv(
  predictor_counts,
  "results/10_rda/RDA_outlier_SNP_counts_by_predictor.csv",
  row.names = FALSE
)

print(axis_counts)
print(predictor_counts)





















############################################################
## Mahalanobis distance histogram
############################################################
#
#mahal_plot <- ggplot(mahal_df, aes(x = Mahalanobis_Dist)) +
#  geom_histogram(
#    binwidth = 1,
#    fill = "steelblue",
#    alpha = 0.6,
#    color = "black"
#  ) +
#  theme_bw() +
#  labs(
#    title = "RDA Mahalanobis Distances",
#    subtitle = "Capblancq/Luu-style method, K = 2, q < 1e-4",
#    x = "Mahalanobis distance",
#    y = "Count"
#  )
#
#if (nrow(outlier_snps) > 0) {
#  mahal_plot <- mahal_plot +
#    geom_vline(
#      xintercept = min(outlier_snps$Mahalanobis_Dist),
#      color = "red",
#      linetype = "dashed"
#    )
#}
#
#print(mahal_plot)
#
#ggsave(
#  "results/10_rda/rda_capblancq_luu_mahalanobis_histogram_K2_q1e-4.pdf",
#  mahal_plot,
#  width = 7,
#  height = 5
#)




























############################################################
## RDA SNP loading plot with strict outliers highlighted
############################################################

plot_df <- snp_loadings
plot_df$is_outlier <- plot_df$SNP %in% outlier_snps$SNP
env_vectors <- as.data.frame(scores(results_full, display = "bp", choices = c(1, 2), scaling = 3))
env_vectors$Variable <- rownames(env_vectors)
snp_radius <- min(max(abs(plot_df$RDA1), na.rm = TRUE),max(abs(plot_df$RDA2), na.rm = TRUE)) * 0.6
env_radius <- max(sqrt(env_vectors$RDA1^2 + env_vectors$RDA2^2),na.rm = TRUE)

arrow_multiplier <- snp_radius / env_radius


rda_outlier_plot <- ggplot() +
  geom_point(data = plot_df,aes(x = RDA1, y = RDA2),color = "gray70",alpha = 0.15,size = 0.4) +
  geom_point(data = subset(plot_df, is_outlier),aes(x = RDA1, y = RDA2),color = "orange2",alpha = 0.95,size = 1.5) +
  geom_segment(data = env_vectors,aes(x = 0,y = 0,xend = RDA1 * arrow_multiplier,yend = RDA2 * arrow_multiplier),arrow = arrow(length = unit(0.18, "cm")),color = "black",linewidth = 0.6) +
  geom_text_repel(data = env_vectors,aes(x = RDA1 * arrow_multiplier,y = RDA2 * arrow_multiplier,label = Variable),color = "black",size = 4,max.overlaps = Inf) +
  coord_equal() +
  theme_bw() +
  labs(x = "RDA1",y = "RDA2") +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))



print(rda_outlier_plot)


ggsave("results/10_rda/RDA_biplot_outlier_snps_highlight.pdf",rda_outlier_plot,width = 8,height = 6)
ggsave("results/10_rda/RDA_biplot_outlier_snps_highlight.png",rda_outlier_plot,width = 8,height = 6,dpi=600)












## Final clean RDA inside-gene table
library(jsonlite)
library(dplyr)
library(readr)
library(stringr)


outliers <- unique_outlier_snps
outliers$scaffold <- outliers$Chromosome
outliers$position <- as.numeric(outliers$Position)



seq_report <- jsonlite::stream_in(file("data/reference_annotation/ncbi_dataset/data/GCF_029207755.1/sequence_report.jsonl"),verbose = FALSE)
name_map <- seq_report %>%
  dplyr::select(scaffold = genbankAccession,annotated_scaffold = refseqAccession)



outliers_mapped <- outliers %>%
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



for (i in seq_len(nrow(outliers_mapped))) {
  snp <- outliers_mapped[i, ]
  
  hits <- genes %>%
    dplyr::filter(
      seqid == snp$annotated_scaffold,
      start <= snp$position,
      end >= snp$position
    )
  
  if (nrow(hits) > 0) {
    hits$SNP <- snp$SNP
    hits$scaffold <- snp$scaffold
    hits$position <- snp$position
    hits$annotated_scaffold <- snp$annotated_scaffold
    inside_hits[[length(inside_hits) + 1]] <- hits
  }
}


inside_genes_clean <- dplyr::bind_rows(inside_hits)



rda_outlier_genes_final <- inside_genes_clean %>%
  dplyr::select(SNP,scaffold,position,annotated_scaffold,gene_id,gene_name,gene_biotype,description,gene_start = start,gene_end = end,strand) %>%
  dplyr::arrange(scaffold, position, gene_name)



write.csv(rda_outlier_genes_final,"results/10_rda/RDA_outlier_GENES.csv",row.names = FALSE)
cat("RDA outlier SNP-gene rows:", nrow(rda_outlier_genes_final), "\n")
cat("Unique RDA outlier SNPs inside genes:", length(unique(rda_outlier_genes_final$SNP)), "\n")
cat("Unique genes containing RDA outlier SNPs:", length(unique(rda_outlier_genes_final$gene_name)), "\n")














































############################################################
## Second RDA using Capblancq/Luu-style outlier SNPs
############################################################

library(raster)
library(terra)
library(sf)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(viridis)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)


# subset genotype matrix to outlier SNPs
genotypes_capblancq_outliers <- genotypes[
  ,
  colnames(genotypes) %in% outlier_snps$SNP,
  drop = FALSE
]


cat("Outlier genotype matrix dimensions:", dim(genotypes_capblancq_outliers), "\n")



# run second RDA
results_capblancq_outliers <- rda_run(
  gen = genotypes_capblancq_outliers,
  env = env_pc_fixed,
  coords = coords_df[, c("Longitude", "Latitude")],
  model = "full",
  correctGEO = FALSE,
  correctPC = FALSE,
  nPC = 2
)



summary(results_capblancq_outliers)
saveRDS(results_capblancq_outliers,"results/10_rda/RDA_adaptive_landscape.rds")
#write.csv(outlier_snps,"results/10_rda/RDA_outlier_snps.csv",row.names = FALSE)










############################################################
## Plot second RDA: samples colored by ecoregion
############################################################

ecoregion_colours <- c(
  "black", "#E69F00", "#56B4E9", "green", "#F0E442",
  "#0072B2", "purple4", "hotpink", "#999999", "blue",
  "orangered", "#4DAF4A", "#984EA3", "red3", "#FFFF33"
)

ecoregion_levels <- levels(factor(metadata$Ecoregion))

ecoregion_palette <- setNames(
  scales::alpha(ecoregion_colours[seq_along(ecoregion_levels)], 0.85),
  ecoregion_levels
)

site_scores_second <- as.data.frame(
  scores(results_capblancq_outliers, display = "sites", choices = c(1, 2), scaling = 3)
)

site_scores_second$SampleID <- rownames(site_scores_second)
site_scores_second$Ecoregion <- factor(metadata$Ecoregion, levels = names(ecoregion_palette))

env_vectors_second <- as.data.frame(
  scores(results_capblancq_outliers, display = "bp", choices = c(1, 2), scaling = 3)
)

env_vectors_second$Variable <- rownames(env_vectors_second)

sample_radius <- min(
  max(abs(site_scores_second$RDA1), na.rm = TRUE),
  max(abs(site_scores_second$RDA2), na.rm = TRUE)
) * 0.75

env_radius <- max(
  sqrt(env_vectors_second$RDA1^2 + env_vectors_second$RDA2^2),
  na.rm = TRUE
)

arrow_multiplier <- sample_radius / env_radius

second_rda_sample_plot <- ggplot() +
  geom_point(
    data = site_scores_second,
    aes(x = RDA1, y = RDA2, color = Ecoregion),
    size = 3,
    alpha = 0.9
  ) +
  geom_segment(
    data = env_vectors_second,
    aes(x = 0, y = 0, xend = RDA1 * arrow_multiplier, yend = RDA2 * arrow_multiplier),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black",
    linewidth = 0.7
  ) +
  geom_text_repel(
    data = env_vectors_second,
    aes(x = RDA1 * arrow_multiplier, y = RDA2 * arrow_multiplier, label = Variable),
    color = "black",
    size = 4,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = ecoregion_palette, drop = FALSE) +
  coord_equal() +
  theme_bw() +
  labs(
    title = "Adaptively Enriched RDA",
    subtitle = "Capblancq/Luu-style outlier SNPs",
    x = "RDA1",
    y = "RDA2",
    color = "Ecoregion"
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.2, alpha = 1)))

print(second_rda_sample_plot)

ggsave(
  "results/10_rda/second_rda_capblancq_luu_outliers_samples_ecoregion_K6_q1e-4.pdf",
  second_rda_sample_plot,
  width = 9,
  height = 7
)













############################################################
## Plot second RDA: SNP loadings with environmental arrows
############################################################

snp_scores_second <- as.data.frame(
  scores(results_capblancq_outliers, display = "species", choices = c(1, 2), scaling = 3)
)

snp_scores_second$SNP <- rownames(snp_scores_second)

env_vectors_second <- as.data.frame(
  scores(results_capblancq_outliers, display = "bp", choices = c(1, 2), scaling = 3)
)

env_vectors_second$Variable <- rownames(env_vectors_second)

snp_radius <- min(
  max(abs(snp_scores_second$RDA1), na.rm = TRUE),
  max(abs(snp_scores_second$RDA2), na.rm = TRUE)
) * 1.7

env_radius <- max(
  sqrt(env_vectors_second$RDA1^2 + env_vectors_second$RDA2^2),
  na.rm = TRUE
)

arrow_multiplier <- snp_radius / env_radius

second_rda_snp_plot <- ggplot() +
  geom_point(
    data = snp_scores_second,
    aes(x = RDA1, y = RDA2),
    color = "gray65",
    alpha = 0.7,
    size = 2
  ) +
  geom_segment(
    data = env_vectors_second,
    aes(
      x = 0,
      y = 0,
      xend = RDA1 * arrow_multiplier,
      yend = RDA2 * arrow_multiplier
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black",
    linewidth = 1
  ) +
  geom_text_repel(
    data = env_vectors_second,
    aes(
      x = RDA1 * arrow_multiplier,
      y = RDA2 * arrow_multiplier,
      label = Variable
    ),
    color = "black",
    size = 4,
    max.overlaps = Inf
  ) +
  theme_minimal() +
  labs(
    x = "RDA1",
    y = "RDA2"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(second_rda_snp_plot)

ggsave(
  "results/10_rda/RDA_biplot_adaptively_enriched_snps.pdf",
  second_rda_snp_plot,
  width = 4.8,
  height = 5
)

ggsave(
  "results/10_rda/RDA_biplot_adaptively_enriched_snps.png",
  second_rda_snp_plot,
  width = 4.8,
  height = 5,
  dpi = 600
)





############################################################
## True adaptive index maps from second RDA
############################################################
############################################################
## Plot adaptive index maps with elevation background
############################################################

library(elevatr)
library(ggnewscale)

# elevation background
study_area <- st_as_sf(
  data.frame(x = c(-130, -115), y = c(30, 47)),
  coords = c("x", "y"),
  crs = 4326
)

elev_raster <- get_elev_raster(
  locations = study_area,
  z = 6,
  clip = "bbox"
)

elev_df <- as.data.frame(elev_raster, xy = TRUE, na.rm = TRUE)
colnames(elev_df) <- c("X", "Y", "Elevation")

elev_min <- min(elev_df$Elevation, na.rm = TRUE)
elev_max <- max(elev_df$Elevation, na.rm = TRUE)
states <- ne_states(country = "United States of America", returnclass = "sf")

bird_data <- st_read("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared/shapefiles/SppDataRequest/SppDataRequest.shp")
bird_data_filtered <- bird_data[bird_data$SCI_NAME == "Chamaea fasciata", ]
bird_data_filtered <- st_transform(bird_data_filtered, crs = 4326)

# environmental loadings from second RDA
rda_env_loadings <- as.data.frame(
  scores(results_capblancq_outliers, display = "bp", choices = c(1, 2))
)

rda_env_loadings <- rda_env_loadings[colnames(env_pc_fixed), ]

rda1_weights <- rda_env_loadings[, "RDA1"]
rda2_weights <- rda_env_loadings[, "RDA2"]

names(rda1_weights) <- rownames(rda_env_loadings)
names(rda2_weights) <- rownames(rda_env_loadings)

# project BIOCLIM rasters into the same climate PC space
names(climate) <- paste0("bio", 1:19)

temp_stack <- climate[[paste0("bio", 1:11)]]
prec_stack <- climate[[paste0("bio", 12:19)]]

temp_pc_rasters <- raster::predict(temp_stack, pca_temp, index = 1:3)
prec_pc_rasters <- raster::predict(prec_stack, pca_prec, index = 1:3)

names(temp_pc_rasters) <- paste0("Temp_PC", 1:3)
names(prec_pc_rasters) <- paste0("Prec_PC", 1:3)

env_pc_rasters <- raster::stack(temp_pc_rasters, prec_pc_rasters)
names(env_pc_rasters) <- colnames(env_pc_fixed)

# standardize raster PCs using sample PC means and SDs
env_centers <- apply(env_pc_fixed, 2, mean, na.rm = TRUE)
env_scales <- apply(env_pc_fixed, 2, sd, na.rm = TRUE)

env_pc_rasters_scaled <- env_pc_rasters

for (v in names(env_pc_rasters_scaled)) {
  env_pc_rasters_scaled[[v]] <- (env_pc_rasters[[v]] - env_centers[v]) / env_scales[v]
}

names(env_pc_rasters_scaled) <- colnames(env_pc_fixed)

# adaptive index = sum(environmental loading * standardized environmental value)
rda1_index_raster <- raster::overlay(
  env_pc_rasters_scaled[[1]],
  env_pc_rasters_scaled[[2]],
  env_pc_rasters_scaled[[3]],
  env_pc_rasters_scaled[[4]],
  env_pc_rasters_scaled[[5]],
  env_pc_rasters_scaled[[6]],
  fun = function(Temp_PC1, Temp_PC2, Temp_PC3, Prec_PC1, Prec_PC2, Prec_PC3) {
    Temp_PC1 * rda1_weights["Temp_PC1"] +
      Temp_PC2 * rda1_weights["Temp_PC2"] +
      Temp_PC3 * rda1_weights["Temp_PC3"] +
      Prec_PC1 * rda1_weights["Prec_PC1"] +
      Prec_PC2 * rda1_weights["Prec_PC2"] +
      Prec_PC3 * rda1_weights["Prec_PC3"]
  }
)

rda2_index_raster <- raster::overlay(
  env_pc_rasters_scaled[[1]],
  env_pc_rasters_scaled[[2]],
  env_pc_rasters_scaled[[3]],
  env_pc_rasters_scaled[[4]],
  env_pc_rasters_scaled[[5]],
  env_pc_rasters_scaled[[6]],
  fun = function(Temp_PC1, Temp_PC2, Temp_PC3, Prec_PC1, Prec_PC2, Prec_PC3) {
    Temp_PC1 * rda2_weights["Temp_PC1"] +
      Temp_PC2 * rda2_weights["Temp_PC2"] +
      Temp_PC3 * rda2_weights["Temp_PC3"] +
      Prec_PC1 * rda2_weights["Prec_PC1"] +
      Prec_PC2 * rda2_weights["Prec_PC2"] +
      Prec_PC3 * rda2_weights["Prec_PC3"]
  }
)

names(rda1_index_raster) <- "RDA1_Index"
names(rda2_index_raster) <- "RDA2_Index"

# mask to Wrentit range
rda1_terra <- terra::rast(rda1_index_raster)
rda2_terra <- terra::rast(rda2_index_raster)

terra::crs(rda1_terra) <- "EPSG:4326"
terra::crs(rda2_terra) <- "EPSG:4326"

bird_range_terra <- terra::vect(bird_data_filtered)
bird_range_terra <- terra::project(bird_range_terra, terra::crs(rda1_terra))

rda1_mask <- terra::mask(terra::crop(rda1_terra, bird_range_terra), bird_range_terra)
rda2_mask <- terra::mask(terra::crop(rda2_terra, bird_range_terra), bird_range_terra)



rda1_df <- as.data.frame(rda1_mask, xy = TRUE, na.rm = TRUE)
rda2_df <- as.data.frame(rda2_mask, xy = TRUE, na.rm = TRUE)

colnames(rda1_df) <- c("X", "Y", "RDA1_Index")
colnames(rda2_df) <- c("X", "Y", "RDA2_Index")

############################################################
## Plot adaptive index maps with elevation background
############################################################

map_rda1 <- ggplot() +
  geom_raster(
    data = elev_df,
    aes(x = X, y = Y, fill = Elevation),
    alpha = 1
  ) +
  scale_fill_gradientn(
    colors = c("white", "white", "gray30"),
    values = scales::rescale(c(elev_min, 500, elev_max)),
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  geom_raster(
    data = rda1_df,
    aes(x = X, y = Y, fill = RDA1_Index),
    alpha = 1
  ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "RDA1 Index"
  ) +
  geom_sf(data = states, fill = NA, color = "black", linewidth = 0.1) +
  coord_sf(xlim = c(-130, -115), ylim = c(30, 47), expand = FALSE) +
  theme_void() +
  theme(legend.position = "right")

map_rda2 <- ggplot() +
  geom_raster(
    data = elev_df,
    aes(x = X, y = Y, fill = Elevation),
    alpha = 1
  ) +
  scale_fill_gradientn(
    colors = c("white", "white", "gray30"),
    values = scales::rescale(c(elev_min, 500, elev_max)),
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  geom_raster(
    data = rda2_df,
    aes(x = X, y = Y, fill = RDA2_Index),
    alpha = 1
  ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "RDA2 Index"
  ) +
  geom_sf(data = states, fill = NA, color = "black", linewidth = 0.1) +
  coord_sf(xlim = c(-130, -115), ylim = c(30, 47), expand = FALSE) +
  theme_void() +
  theme(legend.position = "right")

final_true_adaptive_map <- map_rda1 + map_rda2 + patchwork::plot_layout(ncol = 2)

print(final_true_adaptive_map)

ggsave(
  "results/10_rda/true_adaptive_index_maps_capblancq_luu_outliers_K6_q0.001_RDA1_RDA2_with_elevation.pdf",
  final_true_adaptive_map,
  width = 12,
  height = 6
)

ggsave(
  "results/10_rda/true_adaptive_index_maps_capblancq_luu_outliers_K6_q0.001_RDA1_RDA2_with_elevation.png",
  final_true_adaptive_map,
  width = 12,
  height = 6,
  dpi = 600
)
































############################################################
## Plot adaptive index maps with elevation background
## Long legends
############################################################

library(grid)

map_rda1 <- ggplot() +
  geom_raster(
    data = elev_df,
    aes(x = X, y = Y, fill = Elevation),
    alpha = 1
  ) +
  scale_fill_gradientn(
    colors = c("white", "white", "gray30"),
    values = scales::rescale(c(elev_min, 500, elev_max)),
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  geom_raster(
    data = rda1_df,
    aes(x = X, y = Y, fill = RDA1_Index),
    alpha = 1
  ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "RDA1 Index",
    guide = guide_colorbar(
      barheight = unit(9, "cm"),
      barwidth = unit(0.6, "cm"),
      title.position = "top"
    )
  ) +
  geom_sf(data = states, fill = NA, color = "black", linewidth = 0.1) +
  coord_sf(xlim = c(-130, -115), ylim = c(30, 47), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 14)
  )

map_rda2 <- ggplot() +
  geom_raster(
    data = elev_df,
    aes(x = X, y = Y, fill = Elevation),
    alpha = 1
  ) +
  scale_fill_gradientn(
    colors = c("white", "white", "gray30"),
    values = scales::rescale(c(elev_min, 500, elev_max)),
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  geom_raster(
    data = rda2_df,
    aes(x = X, y = Y, fill = RDA2_Index),
    alpha = 1
  ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "RDA2 Index",
    guide = guide_colorbar(
      barheight = unit(9, "cm"),
      barwidth = unit(0.6, "cm"),
      title.position = "top"
    )
  ) +
  geom_sf(data = states, fill = NA, color = "black", linewidth = 0.1) +
  coord_sf(xlim = c(-130, -115), ylim = c(30, 47), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 14)
  )

final_true_adaptive_map <- map_rda1 + map_rda2 + patchwork::plot_layout(ncol = 2)

print(final_true_adaptive_map)

ggsave(
  "results/10_rda/true_adaptive_index_maps_capblancq_luu_outliers_K6_q0.001_RDA1_RDA2_with_elevation.pdf",
  final_true_adaptive_map,
  width = 12,
  height = 6
)

ggsave(
  "results/10_rda/true_adaptive_index_maps_capblancq_luu_outliers_K6_q0.001_RDA1_RDA2_with_elevation.png",
  final_true_adaptive_map,
  width = 12,
  height = 6,
  dpi = 600
)



















