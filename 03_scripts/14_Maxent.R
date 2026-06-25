#
#climate modeling for black-throated grey warbler
#(c) Michael Anthony Tofflemire - 2022
#install.packages('dismo')
#
##load packages
install.packages("sp")
install.packages("spdep")
install.packages("raster")
install.packages("devtools")
install_github("danlwarren/ENMTools")
install.packages("ENMTools")
install.packages("ENMeval")
install.packages("ecospat")
install.packages("spatstat")
search(ENMTool)
install.packages("pheatmap")
install.packages("BiodiversityR")
install.packages("geosphere")
install.packages("tcltk")
install.packages("ROCR")
packageVersion("ecospat")
library(pheatmap)
library(dismo)
library(raster)
library(sp)
library(rJava)
library(sf)
library(spdep)
library(vegan)
library(devtools)
library(ENMTools)
library(ENMeval)
library(ecospat)
library(ggplot2)
library(tidyverse)
library(spatstat) 
library(terra)
library(dismo)
library(raster)
library(sp)
library(rJava)
library(sf)
library(raster)
library(BiodiversityR)
library(geosphere)
library(tcltk)
library(BiodiversityR)
library(ROCR)
library(ade4)
library(maps)



install.packages("rJava", type = "binary")
library(rJava)
.jinit(); .jcall("java/lang/System","S","getProperty","java.version")
# should print something like "17.x"

#need to start by analyzing data across all Setophaga individuals in dataset
#load climate data
files <- list.files('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared/layers/climate/bio_2-5m_bil',full.names = T)
files <- files[grepl("\\.bil",files)]
climate <- stack(files)



#load occurrence data
#Adjust for spatial autocorrelation 
sp <- read.csv('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/data/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv')
locs <- sp[,c("long","lat")]
locs <- locs[!duplicated(locs), ]
locs <- ensemble.spatialThin(locs, thin.km = 1) #use this to thin based on distance threshold





#crop climate data to make easier to handle 
longbounds <- c(min(sp$long)-6,max(sp$long)+5)
latbounds <- c(min(sp$lat)-7,max(sp$lat)+5)
ext <- extent(c(longbounds,latbounds))
climate <- crop(climate,ext)



# Run intial model without taking out variables
model <- maxent(climate, locs)
pred <- predict(model, climate)
model

################
# Analze results
################
# View variable contributions
variable_contributions <- model@results
head(variable_contributions)


# Access the full results table
model_results <- model@results
print(model_results)


# Plot response curves for each variable
response(model)
# Plot the predicted distribution
plot(pred, main="Predicted Species Distribution")
# Access AUC and other evaluation metrics
auc <- model@evaluation@auc
print(paste("AUC:", auc))


# Save predicted MaxEnt distribution as GeoTIFF
writeRaster(
  pred,
  filename = '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/projects/chamaea_genomics/out/BTGW_maxent_prediction.tif',
  format = "GTiff",
  overwrite = TRUE
)




#Need to now clean up your data based on the preliminary analysis
#check for correlations among variables using Pearson's correlation coefficient
variables <- extract(climate, locs, by = c("lat", "long"))
cor_plot <- cor(variables, method = "pearson")
heatmap(cor_plot, main = "heatmap of Pearson's r values")
# View the correlation matrix directly
# Find pairs of variables with a correlation greater than 0.8
high_cor_pairs <- which(abs(cor_plot) > 0.9, arr.ind = TRUE)


# Filter out self-correlations
high_cor_pairs <- high_cor_pairs[high_cor_pairs[, 1] != high_cor_pairs[, 2], ]

# Display the variable pairs and their correlation values
for (i in 1:nrow(high_cor_pairs)) {
  var1 <- rownames(cor_plot)[high_cor_pairs[i, 1]]
  var2 <- colnames(cor_plot)[high_cor_pairs[i, 2]]
  cor_value <- cor_plot[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
  cat(sprintf("Variables: %s and %s have a correlation of %.2f\n", var1, var2, cor_value))
}

# List of variables you want to remove
remove_vars <- c("bio4", "bio5", "bio8", "bio10", "bio11", "bio12","bio13", "bio16", "bio17")

# Get the names of all layers in the climate stack
layer_names <- names(climate)

# Find the indices of the layers that are not in the remove_vars list
keep_indices <- which(!layer_names %in% remove_vars)

# Subset the climate stack to keep only the layers not in remove_vars
climate_filtered <- climate[[keep_indices]]

# Check the remaining layers in the filtered climate stack
print(names(climate_filtered))







#Process Climate model for contemporary distribution
# Load climate data
files <- list.files("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/bio_2-5m_bil", full.names = TRUE)
files <- files[grepl("\\.bil", files)]
climate <- stack(files)

# List of variables you want to remove
remove_vars <- c("bio4", "bio5", "bio8", "bio10", "bio11", "bio12","bio13", "bio16", "bio17")

# Get the names of all layers in the climate stack
layer_names <- names(climate)

# Find the indices of the layers that are not in the remove_vars list
keep_indices <- which(!layer_names %in% remove_vars)

# Subset the climate stack to keep only the layers not in remove_vars
climate_filtered <- climate[[keep_indices]]

# Load occurrence data and adjust for spatial autocorrelation
sp <- read.csv("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/occurrence/Cfa-metadata-nsamp159copy.csv")
locs <- sp[, c("long", "lat")]
locs <- locs[!duplicated(locs), ]
locs <- ensemble.spatialThin(locs, thin.km = 1) # Thin based on distance threshold

# Crop the filtered climate data to make it easier to handle
longbounds <- c(min(sp$long) - 5, max(sp$long) + 5)
latbounds <- c(min(sp$lat) -2, max(sp$lat) + 2)
ext <- extent(c(longbounds, latbounds))
climate_filtered <- crop(climate_filtered, ext)

# Proceed with further analysis using the filtered and cropped climate stack
# Fit a modern climate model
model <- maxent(climate_filtered, locs)
pred <- predict(model, climate_filtered)
dev.off()
model

plot(pred)
# Load the maps package correctly
library(maps)
# Add country borders
maps::map("world", add = TRUE, col = "black", lwd = 1)
# Add state borders (for USA)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)
# Add Canadian provinces if needed
maps::map("canada", add = TRUE, col = "black", lwd = 0.5)




#pretty plot for contemporary publication
library(paletteer)
pal <- rev(paletteer_c("grDevices::RdYlBu", 100))
plot(pred, 
     col = pal, 
     main = "Current Niche Model (RdYlBu)",
     legend.args = list(text = 'Suitability', side = 4, font = 2, line = 2.5))
maps::map("world", add = TRUE, col = "black", lwd = 1)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)
library(fields)
library(paletteer)
library(maps)
pal <- rev(paletteer_c("grDevices::RdYlBu", 100))
min_val <- minValue(pred)
max_val <- maxValue(pred)
par(mar = c(4, 4, 4, 4), xaxs = "i", yaxs = "i")

# Plot raster
plot(pred,
     col = pal,
     xlim = c(-125, -115),
     ylim = c(ymin(pred), ymax(pred)),
     legend = FALSE,
     axes = FALSE,
     box = FALSE)
box(lwd = 1.7)

maps::map("world", xlim = c(-125, xmax(pred)), ylim = c(ymin(pred), ymax(pred)),
          add = TRUE, col = "black", lwd = 1)
maps::map("state", xlim = c(-125, xmax(pred)), ylim = c(ymin(pred), ymax(pred)),
          add = TRUE, col = "black", lwd = 0.5)

# ✅ Add legend with proper title and clean ticks
image.plot(legend.only = TRUE,
           zlim = c(min_val, max_val),
           col = pal,
           legend.width = 0.3,            # thinner
           legend.shrink = 0.6,           # shorter
           smallplot = c(0.145, 0.185, 0.18, 0.44),  # position (slightly right and down)
           legend.args = list(side = 3, font = 2, line = 1, cex = 0.8),
           legend.axis = list(at = seq(min_val, max_val, length.out = 10),
                              labels = round(seq(min_val, max_val, length.out = 10), 2),
                              cex.axis = 0.7))








###Process LGM files
# List and load all LGM `.tif` files
files <- list.files("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/cclgmbi_2-5m", 
                    pattern = "\\.tif$", full.names = TRUE)

# Stack the files to create a RasterStack
lgm <- stack(files)

# Get the names of all layers in the LGM stack
layer_names_lgm <- names(lgm)

# Extract the numeric suffixes from the layer names (e.g., "cclgmbi10" -> "10")
lgm_suffixes <- as.numeric(gsub("cclgmbi", "", layer_names_lgm))


# List of variable suffixes you want to remove (numeric form to match the LGM naming)
remove_suffixes <- c(4, 5, 8, 10, 11, 12, 13, 16, 17)

# Find the indices of the layers that are not in the remove_suffixes list
lgm_keep_indices <- which(!lgm_suffixes %in% remove_suffixes)

# Subset the LGM stack to keep only the layers not in remove_suffixes
lgm_filtered <- lgm[[lgm_keep_indices]]

# Check the remaining layers in the filtered LGM stack
print(names(lgm_filtered))


# Ensure the number of layers in lgm_filtered matches the number of layers in climate_filtered
if (length(names(lgm_filtered)) == length(names(climate_filtered))) {
  # Assign the names of the filtered climate stack to the LGM stack
  names(lgm_filtered) <- names(climate_filtered)
} else {
  stop("The number of layers in the filtered LGM stack does not match the filtered climate stack.")
}

# Check the renamed layers to confirm
print(names(lgm_filtered))
lgm <- stack(lgm_filtered)
#crop climate data
#I'm cropping it a little differently than before to get the southern extent
#You will need to crop all projections
lgm <- crop(lgm, ext)
lgm_pred <- predict(model,lgm)
plot(lgm_pred)

diff <- lgm_pred - pred
plot(difference)
library(RColorBrewer)
library(viridis)
library(paletteer)
# Create the RdYlBu palette with 30 colors
pal2 <- rev(paletteer_c("grDevices::RdYlBu", 30))
# Example usage in a plot
plot(lgm_pred, 
     col = pal2, 
     main = "LGM Projection: RdYlBu Color Ramp",
     legend.args = list(text = 'Probability', side = 4, font = 2, line = 2.5))
# Load the maps package correctly
library(maps)
# Add country borders
maps::map("world", add = TRUE, col = "black", lwd = 1)
# Add state borders (for USA)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)
# Add Canadian provinces if needed
maps::map("canada", add = TRUE, col = "black", lwd = 0.5)
library(fields)
library(paletteer)
library(maps)
pal <- rev(paletteer_c("grDevices::RdYlBu", 100))
min_val <- minValue(lgm_pred)
max_val <- maxValue(lgm_pred)
par(mar = c(4, 4, 4, 4), xaxs = "i", yaxs = "i")

# Plot raster
plot(lgm_pred,
     col = pal,
     xlim = c(-125, -115),
     ylim = c(ymin(lgm_pred), ymax(lgm_pred)),
     legend = FALSE,
     axes = FALSE,
     box = FALSE)
box(lwd = 1.7)

maps::map("world", xlim = c(-125, xmax(pred)), ylim = c(ymin(lgm_pred), ymax(lgm_pred)),
          add = TRUE, col = "black", lwd = 1)
maps::map("state", xlim = c(-125, xmax(pred)), ylim = c(ymin(lgm_pred), ymax(lgm_pred)),
          add = TRUE, col = "black", lwd = 0.5)

# ✅ Add legend with proper title and clean ticks
image.plot(legend.only = TRUE,
           zlim = c(min_val, max_val),
           col = pal,
           legend.width = 0.3,            # thinner
           legend.shrink = 0.6,           # shorter
           smallplot = c(0.145, 0.185, 0.18, 0.44),  # position (slightly right and down)
           legend.args = list(side = 3, font = 2, line = 1, cex = 0.8),
           legend.axis = list(at = seq(min_val, max_val, length.out = 10),
                              labels = round(seq(min_val, max_val, length.out = 10), 2),
                              cex.axis = 0.7))





###############################






















# List and load all LGM `.tif` files
files <- list.files("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/lig_30s_bio", 
                    pattern = "\\.bil$", full.names = TRUE)

# Stack the files to create a RasterStack
lig <- stack(files)

# Get the names of all layers in the LGM stack
layer_names_lig <- names(lig)

# Extract the numeric suffixes from the layer names (e.g., "cclgmbi10" -> "10")
lig_suffixes <- as.numeric(gsub("lig_30s_bio_", "", layer_names_lig))

# List of variable suffixes you want to remove (numeric form to match the LGM naming)
remove_suffixes <- c(4, 5, 8, 10, 11, 12, 13, 16, 17)

# Find the indices of the layers that are not in the remove_suffixes list
lig_keep_indices <- which(!lig_suffixes %in% remove_suffixes)

# Subset the LGM stack to keep only the layers not in remove_suffixes
lig_filtered <- lig[[lig_keep_indices]]

# Check the remaining layers in the filtered LGM stack
print(names(lig_filtered))
lig_filtered

# Ensure the number of layers in lgm_filtered matches the number of layers in climate_filtered
if (length(names(lig_filtered)) == length(names(climate_filtered))) {
  # Assign the names of the filtered climate stack to the LGM stack
  names(lig_filtered) <- names(climate_filtered)
} else {
  stop("The number of layers in the filtered LGM stack does not match the filtered climate stack.")
}

# Check the renamed layers to confirm
print(names(lig_filtered))


#####
print(lig_suffixes)
######
lig <- stack(lig_filtered)


#crop climate data
#I'm cropping it a little differently than before to get the southern extent
#You will need to crop all projections
lig <- crop(lig, ext)
lig_pred <- predict(model,lig)
plot(lig_pred)
library(fields)
library(paletteer)
library(maps)
pal <- rev(paletteer_c("grDevices::RdYlBu", 100))
min_val <- minValue(lig_pred)
max_val <- maxValue(lig_pred)
par(mar = c(4, 4, 4, 4), xaxs = "i", yaxs = "i")

# Plot raster
plot(lig_pred,
     col = pal,
     xlim = c(-125, -115),
     ylim = c(ymin(lig_pred), ymax(lig_pred)),
     legend = FALSE,
     axes = FALSE,
     box = FALSE)
box(lwd = 1.7)

maps::map("world", xlim = c(-125, xmax(pred)), ylim = c(ymin(lig_pred), ymax(lig_pred)),
          add = TRUE, col = "black", lwd = 1)
maps::map("state", xlim = c(-125, xmax(pred)), ylim = c(ymin(lig_pred), ymax(lig_pred)),
          add = TRUE, col = "black", lwd = 0.5)

# ✅ Add legend with proper title and clean ticks
image.plot(legend.only = TRUE,
           zlim = c(min_val, max_val),
           col = pal,
           legend.width = 0.3,            # thinner
           legend.shrink = 0.6,           # shorter
           smallplot = c(0.145, 0.185, 0.18, 0.44),  # position (slightly right and down)
           legend.args = list(side = 3, font = 2, line = 1, cex = 0.8),
           legend.axis = list(at = seq(min_val, max_val, length.out = 10),
                              labels = round(seq(min_val, max_val, length.out = 10), 2),
                              cex.axis = 0.7))





####
# List and load all LGM `.tif` files
files <- list.files("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/mrmidbi_2-5m", 
                    pattern = "\\.tif$", full.names = TRUE)
# Stack the files to create a RasterStack
mh <- stack(files)
# Get the names of all layers in the LGM stack
layer_names_mh <- names(mh)
# Extract the numeric suffixes from the layer names (e.g., "cclgmbi10" -> "10")
mh_suffixes <- as.numeric(gsub("mrmidbi", "", layer_names_mh))


# List of variable suffixes you want to remove (numeric form to match the LGM naming)
remove_suffixes <- c(4, 5, 8, 10, 11, 12, 13, 16, 17)

# Find the indices of the layers that are not in the remove_suffixes list
mh_keep_indices <- which(!mh_suffixes %in% remove_suffixes)

# Subset the LGM stack to keep only the layers not in remove_suffixes
mh_filtered <- mh[[mh_keep_indices]]

# Check the remaining layers in the filtered LGM stack
print(names(mh_filtered))

# Ensure the number of layers in lgm_filtered matches the number of layers in climate_filtered
if (length(names(mh_filtered)) == length(names(climate_filtered))) {
  # Assign the names of the filtered climate stack to the LGM stack
  names(mh_filtered) <- names(climate_filtered)
} else {
  stop("The number of layers in the filtered LGM stack does not match the filtered climate stack.")
}

# Check the renamed layers to confirm
print(names(mh_filtered))
mh <- stack(mh_filtered)
#crop climate data
#I'm cropping it a little differently than before to get the southern extent
#You will need to crop all projections
mh <- crop(mh, ext)
mh_pred <- predict(model,mh)
plot(mh_pred)
library(fields)
library(paletteer)
library(maps)
pal <- rev(paletteer_c("grDevices::RdYlBu", 100))
min_val <- minValue(mh_pred)
max_val <- maxValue(mh_pred)
par(mar = c(4, 4, 4, 4), xaxs = "i", yaxs = "i")

# Plot raster
plot(mh_pred,
     col = pal,
     xlim = c(-125, -115),
     ylim = c(ymin(mh_pred), ymax(mh_pred)),
     legend = FALSE,
     axes = FALSE,
     box = FALSE)
box(lwd = 1.7)

maps::map("world", xlim = c(-125, xmax(pred)), ylim = c(ymin(mh_pred), ymax(mh_pred)),
          add = TRUE, col = "black", lwd = 1)
maps::map("state", xlim = c(-125, xmax(pred)), ylim = c(ymin(mh_pred), ymax(mh_pred)),
          add = TRUE, col = "black", lwd = 0.5)

# ✅ Add legend with proper title and clean ticks
image.plot(legend.only = TRUE,
           zlim = c(min_val, max_val),
           col = pal,
           legend.width = 0.3,            # thinner
           legend.shrink = 0.6,           # shorter
           smallplot = c(0.145, 0.185, 0.18, 0.44),  # position (slightly right and down)
           legend.args = list(side = 3, font = 2, line = 1, cex = 0.8),
           legend.axis = list(at = seq(min_val, max_val, length.out = 10),
                              labels = round(seq(min_val, max_val, length.out = 10), 2),
                              cex.axis = 0.7))





###############################





















































# This code makes separate ENMs for each structured population
# In this case, it's for "North" and "South" populations

# Load necessary libraries
library(dismo)
library(raster)
library(sp)
library(ENMeval)

# Load climate data
files <- list.files("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Chamaea_genomics/data/climate/bio_2-5m_bil", full.names = TRUE)
files <- files[grepl("\\.bil", files)]
climate <- stack(files)

# Load coordinates and population assignments
# Load coordinates with sample ID
coords <- read.table(
  "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Chamaea_genomics/data/41-Chamaea.coords.txt",
  header = FALSE,
  col.names = c("sample_id", "long", "lat")
)

# Check the first few rows to confirm successful loading
head(coords)
# Load population assignments with sample IDs
populations <- read.table(
  "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Chamaea_genomics/data/41-Chamaea.pop.reordered.txt",
  header = FALSE,
  col.names = c("sample_id", "population")
)

# Check the first few rows to confirm successful loading
head(populations)

# Combine coordinates with population data
data <- cbind(coords, populations)

# Split into two populations
pop1 <- subset(data, population == "North")
pop2 <- subset(data, population == "South")

library(sp)

# Define a function for MaxEnt modeling with consistent extent
run_maxent <- function(locs, climate_data, crop_extent, output_name) {
  # Remove duplicate locations (optional, but keeps the data clean)
  locs <- locs[!duplicated(locs), ]
  
  # Crop climate data using the pre-defined extent
  cropped_climate <- crop(climate_data, crop_extent)
  
  # Convert locations to SpatialPoints
  locs_sp <- SpatialPoints(locs[, c("long", "lat")], proj4string = crs(cropped_climate))
  
  # Run MaxEnt model
  model <- maxent(cropped_climate, locs_sp)
  pred <- predict(model, cropped_climate)
  
  # Save predictions
  writeRaster(pred, paste0(output_name, "_prediction.tif"), format = "GTiff", overwrite = TRUE)
  
  return(list(model = model, prediction = pred))
}


# Combine all locations
all_locs <- rbind(pop1[, c("long", "lat")], pop2[, c("long", "lat")])

# Calculate overall bounds for cropping
longbounds <- c(min(all_locs$long) - 5, max(all_locs$long) + 5)
latbounds <- c(min(all_locs$lat) - 2, max(all_locs$lat) + 2)
crop_extent <- extent(c(longbounds, latbounds))

# Run MaxEnt for Population 1
result_pop1 <- run_maxent(pop1, climate, crop_extent, "pop1")

# Run MaxEnt for Population 2
result_pop2 <- run_maxent(pop2, climate, crop_extent, "pop2")


# Output summary
result_pop1$model
result_pop2$model


# Load necessary libraries
library(raster)
library(ggplot2)
library(sf)

# Load prediction rasters
pop1_pred <- raster("pop1_prediction.tif")
pop2_pred <- raster("pop2_prediction.tif")
plot(pop1_pred)
plot(pop2_pred)

dev.off()
diff <- pop1_pred-pop2_pred
plot(diff)






















# Load climate data
files <- list.files("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/bio_2-5m_bil", full.names = TRUE)
files <- files[grepl("\\.bil", files)]
climate <- stack(files)

# Remove variables you don't want (must match what's removed from LIG)
remove_vars <- c("bio4", "bio5", "bio7", "bio8", "bio10", "bio11", "bio12", "bio13", "bio16", "bio17", "bio19")
climate_filtered <- climate[[which(!names(climate) %in% remove_vars)]]

# Load occurrence data and thin
sp <- read.csv("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/occurrence/Cfa-metadata-nsamp159copy.csv")
locs <- sp[, c("long", "lat")]
locs <- locs[!duplicated(locs), ]
locs <- ensemble.spatialThin(locs, thin.km = 1)

# Crop to relevant extent
longbounds <- c(min(sp$long) - 5, max(sp$long) + 5)
latbounds <- c(min(sp$lat) - 2, max(sp$lat) + 2)
ext <- extent(c(longbounds, latbounds))
climate_filtered <- crop(climate_filtered, ext)

# Fit the MaxEnt model
model <- maxent(climate_filtered, locs)

# (Optional but recommended) Save the model to disk
saveRDS(model, "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/out/climate-models/model_trained_filtered_vars.rds")
model <- readRDS("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/out/climate-models/model_trained_filtered_vars.rds")

library(raster)

library(raster)

# 1. Load LIG layers
lig_dir <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/lig_30s_bio"
files <- list.files(lig_dir, pattern = "\\.bil$", full.names = TRUE)
lig <- stack(files)

# 2. Extract numeric suffixes from LIG names
lig_suffixes <- as.numeric(gsub("lig_30s_bio_", "", names(lig)))

# 3. Remove unwanted variables
remove_suffixes <- c(4, 5, 7, 8, 10, 11, 12, 13, 16, 17, 19)
keep_indices <- which(!lig_suffixes %in% remove_suffixes)
lig_filtered <- lig[[keep_indices]]

# 4. Match names with climate_filtered (assuming you already filtered it with these same vars removed)
names(lig_filtered) <- names(climate_filtered)

# Crop LIG layers to same extent as modern data
lig_filtered <- crop(lig_filtered, ext)
# 5. Project model
lig_pred <- predict(model, lig_filtered)

# 6. Optional: Plot
plot(lig_pred, main = "Projected Suitability – LIG")

# 7. Optional: Save to file
writeRaster(lig_pred,
            filename = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/LIG_maxent_projection.tif",
            format = "GTiff", overwrite = TRUE
)













library(raster)
library(dismo)

# === 1. Load modern climate data ===
mod_files <- list.files(
  "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/bio_2-5m_bil",
  pattern = "\\.bil$",
  full.names = TRUE
)
climate <- stack(mod_files)

# === 2. Remove selected variables ===
remove_vars <- c("bio4", "bio5", "bio8", "bio10", "bio11", "bio12", "bio13", "bio16", "bio17")
climate_filtered <- climate[[which(!names(climate) %in% remove_vars)]]

# === 3. Load and thin occurrence data ===
sp <- read.csv("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/occurrence/Cfa-metadata-nsamp159copy.csv")
locs <- sp[, c("long", "lat")]
locs <- locs[!duplicated(locs), ]
locs <- ensemble.spatialThin(locs, thin.km = 1)

# === 4. Crop to relevant extent ===
longbounds <- c(min(sp$long) - 5, max(sp$long) + 5)
latbounds <- c(min(sp$lat) - 2, max(sp$lat) + 2)
ext <- extent(c(longbounds, latbounds))
climate_filtered <- crop(climate_filtered, ext)

# === 5. Train MaxEnt model ===
model <- maxent(climate_filtered, locs)

# (Optional) Save model
saveRDS(model, "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/out/climate-models/model_filtered_vars.rds")

# === 6. Load and process LIG layers ===
lig_dir <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/lig_30s_bio"
lig_files <- list.files(lig_dir, pattern = "\\.bil$", full.names = TRUE)
lig <- stack(lig_files)

# === 7. Extract LIG variable numbers and remove matching ones ===
lig_suffixes <- as.numeric(gsub("lig_30s_bio_", "", names(lig)))
remove_suffixes <- as.numeric(gsub("bio", "", remove_vars))
keep_indices <- which(!lig_suffixes %in% remove_suffixes)
lig_filtered <- lig[[keep_indices]]

# === 8. Match LIG layer names to climate_filtered ===
names(lig_filtered) <- names(climate_filtered)

# === 9. Crop LIG to same extent ===
lig_filtered <- crop(lig_filtered, ext)

# === 10. Project model onto LIG data ===
lig_pred <- predict(model, lig_filtered)

# === 11. Plot and save ===
plot(lig_pred, main = "Projected Suitability – LIG (Filtered Bioclim Vars)")

writeRaster(
  lig_pred,
  filename = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/LIG_maxent_projection_filtered_vars.tif",
  format = "GTiff",
  overwrite = TRUE
)






















library(raster)
library(dismo)
library(viridis)
library(RColorBrewer)

### === 1. Load and Filter Modern Climate ===
mod_files <- list.files(
  "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/bio_2-5m_bil",
  pattern = "\\.bil$", full.names = TRUE
)
climate <- stack(mod_files)

remove_vars <- c("bio4", "bio5", "bio8", "bio10", "bio11", "bio12", "bio13", "bio16", "bio17")
keep_indices <- which(!names(climate) %in% remove_vars)
climate_filtered <- climate[[keep_indices]]

### === 2. Load and Thin Occurrence Data ===
sp <- read.csv("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/occurrence/Cfa-metadata-nsamp159copy.csv")
locs <- sp[, c("long", "lat")]
locs <- locs[!duplicated(locs), ]
locs <- ensemble.spatialThin(locs, thin.km = 1)

### === 3. Crop Modern Climate ===
longbounds <- c(min(sp$long) - 5, max(sp$long) + 5)
latbounds <- c(min(sp$lat) - 2, max(sp$lat) + 2)
ext <- extent(c(longbounds, latbounds))
climate_filtered <- crop(climate_filtered, ext)

### === 4. Train MaxEnt Model ===
model <- maxent(climate_filtered, locs)
current_pred <- predict(model, climate_filtered)

### === 5. Load and Process LGM ===
lgm_files <- list.files(
  "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/cclgmbi_2-5m",
  pattern = "\\.tif$", full.names = TRUE
)
lgm <- stack(lgm_files)
lgm_suffixes <- as.numeric(gsub("cclgmbi", "", names(lgm)))
remove_suffixes <- c(4, 5, 8, 10, 11, 12, 13, 16, 17)
lgm_keep <- which(!lgm_suffixes %in% remove_suffixes)
lgm_filtered <- lgm[[lgm_keep]]
names(lgm_filtered) <- names(climate_filtered)
lgm_filtered <- crop(lgm_filtered, ext)
lgm_pred <- predict(model, lgm_filtered)

### === 6. Load and Process LIG ===
lig_files <- list.files(
  "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/lig_30s_bio",
  pattern = "\\.bil$", full.names = TRUE
)
lig <- stack(lig_files)
lig_suffixes <- as.numeric(gsub("lig_30s_bio_", "", names(lig)))
lig_keep <- which(!lig_suffixes %in% remove_suffixes)
lig_filtered <- lig[[lig_keep]]
names(lig_filtered) <- names(climate_filtered)
lig_filtered <- crop(lig_filtered, ext)
lig_pred <- predict(model, lig_filtered)

### === 7. Plot Projections ===
par(mfrow = c(3, 1), mar = c(3, 3, 2, 1))

plot(current_pred, main = "Current Suitability")
plot(lgm_pred, main = "LGM Suitability")
plot(lig_pred, main = "LIG Suitability")

### === 8. Plot Differences ===
par(mfrow = c(1, 2))

plot(lgm_pred - current_pred,
     main = "LGM - Current",
     col = rev(brewer.pal(11, "RdBu")),
     zlim = c(-1, 1))

plot(lig_pred - current_pred,
     main = "LIG - Current",
     col = rev(brewer.pal(11, "RdBu")),
     zlim = c(-1, 1))

### === 9. Save Rasters ===
writeRaster(current_pred,
            "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/current_maxent_projection.tif",
            format = "GTiff", overwrite = TRUE)

writeRaster(lgm_pred,
            "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/LGM_maxent_projection.tif",
            format = "GTiff", overwrite = TRUE)

writeRaster(lig_pred,
            "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/LIG_maxent_projection.tif",
            format = "GTiff", overwrite = TRUE)

writeRaster(lgm_pred - current_pred,
            "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/diff_LGM_minus_current.tif",
            format = "GTiff", overwrite = TRUE)

writeRaster(lig_pred - current_pred,
            "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/diff_LIG_minus_current.tif",
            format = "GTiff", overwrite = TRUE)











library(raster)
library(grid)
library(gridExtra)
library(paletteer)
library(maps)
library(png)

# Custom function to generate a grob from raster
plot_raster_to_grob <- function(r, col_palette) {
  tmpfile <- tempfile(fileext = ".png")
  png(tmpfile, width = 800, height = 800, res = 150)
  par(mar = c(0, 0, 0, 0))
  plot(r, col = col_palette, axes = TRUE, box = TRUE, legend = TRUE)
  maps::map("world", add = TRUE, col = "black", lwd = 0.7)
  maps::map("state", add = TRUE, col = "black", lwd = 0.4)
  dev.off()
  img <- readPNG(tmpfile)
  file.remove(tmpfile)
  rasterGrob(img, interpolate = TRUE)
}

# Color palette
pal2 <- rev(paletteer_c("grDevices::RdYlBu", 30))

# Create grobs
g1 <- plot_raster_to_grob(current_pred, pal2)
g2 <- plot_raster_to_grob(lgm_pred, pal2)
g3 <- plot_raster_to_grob(lig_pred, pal2)

# Combine them side by side
grid.newpage()
grid.arrange(g1, g2, g3, ncol = 3, widths = unit(rep(1, 3), "null"))
























# Load necessary libraries
library(dismo)
library(raster)
library(sp)
library(maps)

# ---- Load and filter climate data ----
files <- list.files("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/bio_2-5m_bil", 
                    pattern = "\\.bil$", full.names = TRUE)
climate <- stack(files)

# Drop selected variables (same ones as earlier)
remove_vars <- c("bio4", "bio5", "bio8", "bio10", "bio11", "bio12","bio13", "bio16", "bio17")
layer_names <- names(climate)
keep_indices <- which(!layer_names %in% remove_vars)
climate_filtered <- climate[[keep_indices]]

# ---- Load and merge sample coordinates and population data ----
coords <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/ccgp/qc/41-Chamaea.coords.txt",
                     header = FALSE, col.names = c("sample_id", "long", "lat"))

populations <- read.table("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/gvcfs/arch/41-Chamaea.pop.reordered_OR.txt",
                          header = FALSE, col.names = c("sample_id", "population"))

data <- merge(coords, populations, by = "sample_id")

# ---- Subset North and South ----
north_locs <- subset(data, population == "North")[, c("long", "lat")]
south_locs <- subset(data, population == "South")[, c("long", "lat")]

# ---- Crop extent based on all points ----
all_locs <- rbind(north_locs, south_locs)
longbounds <- c(min(all_locs$long) - 5, max(all_locs$long) + 5)
latbounds <- c(min(all_locs$lat) - 2, max(all_locs$lat) + 2)
ext <- extent(c(longbounds, latbounds))
climate_filtered <- crop(climate_filtered, ext)

# ---- Format locations as SpatialPoints ----
north_points <- SpatialPoints(north_locs, proj4string = crs(climate_filtered))
south_points <- SpatialPoints(south_locs, proj4string = crs(climate_filtered))

# ---- Run MaxEnt separately ----
model_north <- maxent(climate_filtered, north_points)
model_south <- maxent(climate_filtered, south_points)

# ---- Predict habitat suitability ----
pred_north <- predict(model_north, climate_filtered)
pred_south <- predict(model_south, climate_filtered)

# ---- Save prediction rasters ----
writeRaster(pred_north, "maxent_north_prediction.tif", format = "GTiff", overwrite = TRUE)
writeRaster(pred_south, "maxent_south_prediction.tif", format = "GTiff", overwrite = TRUE)

# ---- Visualize predictions ----
par(mfrow = c(1, 3), mar = c(0.5, 0.5, 1.5, 0.5), oma = c(0, 0, 0, 0))

pal <- rev(paletteer::paletteer_c("grDevices::RdYlBu", 30))

plot(pred_north, 
     col = pal, 
     main = "North Cluster", 
     axes = FALSE, box = FALSE, legend = FALSE, zlim = c(0,1))
maps::map("world", add = TRUE, col = "black", lwd = 1)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)

plot(pred_south, 
     col = pal, 
     main = "South Cluster", 
     axes = FALSE, box = FALSE, legend = FALSE, zlim = c(0,1))
maps::map("world", add = TRUE, col = "black", lwd = 1)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)

plot(pred_north - pred_south,
     col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
     main = "North - South Difference", 
     axes = FALSE, box = FALSE, legend = FALSE, zlim = c(-1, 1))
maps::map("world", add = TRUE, col = "black", lwd = 1)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)

par(mfrow = c(1,1))







# Load and filter LGM
lgm_files <- list.files(
  "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/cclgmbi_2-5m",
  pattern = "\\.tif$", full.names = TRUE
)
lgm <- stack(lgm_files)

# Match LGM layers to modern filtered layers
lgm_suffixes <- as.numeric(gsub("cclgmbi", "", names(lgm)))
remove_suffixes <- c(4, 5, 8, 10, 11, 12, 13, 16, 17)
keep_lgm <- which(!lgm_suffixes %in% remove_suffixes)
lgm_filtered <- lgm[[keep_lgm]]
names(lgm_filtered) <- names(climate_filtered)
lgm_filtered <- crop(lgm_filtered, ext)






# Load and filter LIG
lig_files <- list.files(
  "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/lig_30s_bio",
  pattern = "\\.bil$", full.names = TRUE
)
lig <- stack(lig_files)

# Match LIG layers to modern filtered layers
lig_suffixes <- as.numeric(gsub("lig_30s_bio_", "", names(lig)))
keep_lig <- which(!lig_suffixes %in% remove_suffixes)
lig_filtered <- lig[[keep_lig]]
names(lig_filtered) <- names(climate_filtered)
lig_filtered <- crop(lig_filtered, ext)







# LGM predictions
lgm_pred_north <- predict(model_north, lgm_filtered)
lgm_pred_south <- predict(model_south, lgm_filtered)

# LIG predictions
lig_pred_north <- predict(model_north, lig_filtered)
lig_pred_south <- predict(model_south, lig_filtered)





writeRaster(lgm_pred_north, "maxent_LGM_north.tif", overwrite=TRUE)
writeRaster(lgm_pred_south, "maxent_LGM_south.tif", overwrite=TRUE)
writeRaster(lig_pred_north, "maxent_LIG_north.tif", overwrite=TRUE)
writeRaster(lig_pred_south, "maxent_LIG_south.tif", overwrite=TRUE)







par(mfrow = c(2, 2), mar = c(0.5, 0.5, 1.5, 0.5))

# LGM - North
plot(lgm_pred_north, 
     col = pal, 
     main = "LGM - North", 
     axes = FALSE, box = FALSE, legend = FALSE, zlim = c(0, 1))
maps::map("world", add = TRUE, col = "black", lwd = 1)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)

# LGM - South
plot(lgm_pred_south, 
     col = pal, 
     main = "LGM - South", 
     axes = FALSE, box = FALSE, legend = FALSE, zlim = c(0, 1))
maps::map("world", add = TRUE, col = "black", lwd = 1)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)

# LIG - North
plot(lig_pred_north, 
     col = pal, 
     main = "LIG - North", 
     axes = FALSE, box = FALSE, legend = FALSE, zlim = c(0, 1))
maps::map("world", add = TRUE, col = "black", lwd = 1)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)

# LIG - South
plot(lig_pred_south, 
     col = pal, 
     main = "LIG - South", 
     axes = FALSE, box = FALSE, legend = FALSE, zlim = c(0, 1))
maps::map("world", add = TRUE, col = "black", lwd = 1)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)

par(mfrow = c(1,1))
















































# List and load all LGM `.tif` files
files <- list.files("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/climate/mrmidbi_2-5m", 
                    pattern = "\\.tif$", full.names = TRUE)

# Stack the files to create a RasterStack
mh <- stack(files)

# Get the names of all layers in the LGM stack
layer_names_mh <- names(mh)

# Extract the numeric suffixes from the layer names (e.g., "cclgmbi10" -> "10")
mh_suffixes <- as.numeric(gsub("mrmidbi", "", layer_names_mh))


# List of variable suffixes you want to remove (numeric form to match the LGM naming)
remove_suffixes <- c(4, 5, 8, 10, 11, 12, 13, 16, 17)

# Find the indices of the layers that are not in the remove_suffixes list
mh_keep_indices <- which(!mh_suffixes %in% remove_suffixes)

# Subset the LGM stack to keep only the layers not in remove_suffixes
mh_filtered <- mh[[mh_keep_indices]]

# Check the remaining layers in the filtered LGM stack
print(names(mh_filtered))


# Ensure the number of layers in lgm_filtered matches the number of layers in climate_filtered
if (length(names(mh_filtered)) == length(names(climate_filtered))) {
  # Assign the names of the filtered climate stack to the LGM stack
  names(mh_filtered) <- names(climate_filtered)
} else {
  stop("The number of layers in the filtered LGM stack does not match the filtered climate stack.")
}

# Check the renamed layers to confirm
print(names(mh_filtered))







mh <- stack(mh_filtered)


#crop climate data
#I'm cropping it a little differently than before to get the southern extent
#You will need to crop all projections
mh <- crop(mh, ext)
mh_pred <- predict(model,mh)
plot(mh_pred)

diff <- lgm_pred - pred
plot(difference)
library(RColorBrewer)
library(viridis)

# Use viridis color scale
plot(diff, 
     col = viridis(100), 
     main = "Difference Between LGM and Current Prediction",
     xlab = "Longitude",
     ylab = "Latitude")
plot(diff, 
     col = viridis(100, option = "magma"), 
     main = "Difference Between LGM and Current Prediction",
     xlab = "Longitude",
     ylab = "Latitude")




# Create the custom color ramp using adjusted, richer RGB values
richer_colors <- colorRampPalette(c(rgb(60/255, 100/255, 160/255), 
                                    rgb(245/255, 210/255, 100/255), 
                                    rgb(220/255, 40/255, 40/255)))

# Plot the projection with the custom, richer color ramp
plot(mh_pred, 
     col = richer_colors(100), 
     main = "LGM Projection: Richer Color Ramp (Blue to Yellow to Red)",
     legend.args = list(text = 'Probability', side = 4, font = 2, line = 2.5))



library(paletteer)

# Create the RdYlBu palette with 30 colors
pal2 <- rev(paletteer_c("grDevices::RdYlBu", 30))

# Example usage in a plot
plot(mh_pred, 
     col = pal2, 
     main = "LGM Projection: RdYlBu Color Ramp",
     legend.args = list(text = 'Probability', side = 4, font = 2, line = 2.5))








# Load the maps package correctly
library(maps)

# Plot the LGM projection with the custom, richer color ramp
plot(mh_pred, 
     col = richer_colors(100), 
     main = "LGM Projection: Richer Color Ramp (Blue to Yellow to Red)",
     legend.args = list(text = 'Probability', side = 4, font = 2, line = 2.5))

# Add country borders
maps::map("world", add = TRUE, col = "black", lwd = 1)

# Add state borders (for USA)
maps::map("state", add = TRUE, col = "black", lwd = 0.5)

# Add Canadian provinces if needed
maps::map("canada", add = TRUE, col = "black", lwd = 0.5)








library(fields)
library(paletteer)
library(maps)

pal <- rev(paletteer_c("grDevices::RdYlBu", 100))
min_val <- minValue(pred)
max_val <- maxValue(pred)

par(mar = c(4, 4, 4, 4), xaxs = "i", yaxs = "i")

# Plot raster
plot(pred,
     col = pal,
     xlim = c(-125, -115),
     ylim = c(ymin(pred), ymax(pred)),
     legend = FALSE,
     axes = FALSE,
     box = FALSE)
box(lwd = 1.7)

maps::map("world", xlim = c(-125, xmax(pred)), ylim = c(ymin(pred), ymax(pred)),
          add = TRUE, col = "black", lwd = 1)
maps::map("state", xlim = c(-125, xmax(pred)), ylim = c(ymin(pred), ymax(pred)),
          add = TRUE, col = "black", lwd = 0.5)

# ✅ Add legend with proper title and clean ticks
image.plot(legend.only = TRUE,
           zlim = c(min_val, max_val),
           col = pal,
           legend.width = 0.3,            # thinner
           legend.shrink = 0.6,           # shorter
           smallplot = c(0.145, 0.185, 0.18, 0.44),  # position (slightly right and down)
           legend.args = list(side = 3, font = 2, line = 1, cex = 0.8),
           legend.axis = list(at = seq(min_val, max_val, length.out = 10),
                              labels = round(seq(min_val, max_val, length.out = 10), 2),
                              cex.axis = 0.7))


