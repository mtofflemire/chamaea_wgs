## Part 1: Install rEEMSplots
## Check that the current directory contains the rEEMSplots source directory
if (file.exists("/Users/michaeltofflemire/Library/CloudStorage/Dropbox/systems/eems/plotting/rEEMSplots")) {
  install.packages("rEEMSplots", repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}


## Possibly change the working directory with setwd()


# Part 2: Generate graphics
library(rEEMSplots)

# Define the paths
mcmcpath <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/chamaea-eems-analysis/chamaea_eems_autosomal_run1"
plotpath <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/chamaea-eems-analysis/chamaea_eems_autosomal_run1"

# Run the plotting function
eems.plots(mcmcpath, plotpath, longlat = TRUE)

eems.

eems.plots(
  mcmcpath, 
  plotpath, 
  longlat = TRUE,  # Keep lat-long projection
  add.outline = TRUE,  # Add an outline around EEMS grids
  add.grid = TRUE,  # Remove the hexagonal grid for a cleaner map
)


library(rEEMSplots)
library(rEEMSplots)
library(maps)
library(mapdata)
library(rworldmap)
library(rworldxtra)

projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"



eems.plots(
  mcmcpath = mcmcpath,
  plotpath = paste0(plotpath, "-geographic-map"),
  longlat = TRUE,
  projection.in = projection_none,
  projection.out = projection_mercator,
  add.map = TRUE,  # This adds a high-res map
  col.map = "black",  # Change map outline color
  lwd.map = 2,  # Line width for the map
  add.demes=TRUE
)


eems.p


map_world <- getMap(resolution = "high")  # Load high-res world map
map_CA_OR <- map_world[which(map_world@data$SOVEREIGNT == "United States"), ]

eems.plots(
  mcmcpath = mcmcpath,
  plotpath = paste0(plotpath, "-shapefile"),
  longlat = TRUE,
  m.plot.xy = {
    plot(map_CA_OR, col = NA, border = "black", add = TRUE)
  },
  q.plot.xy = {
    plot(map_CA_OR, col = NA, border = "black", add = TRUE)
  }
)






\library(rEEMSplots)
library(maps)
library(mapdata)

# Define paths
mcmcpath <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/chamaea-eems-analysis/chamaea_eems_run3"
plotpath <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/chamaea-eems-analysis/plots"

# Generate EEMS plots with California & Oregon state outlines
eems.plots(
  mcmcpath = mcmcpath,
  plotpath = paste0(plotpath, "-state-outlines"),
  longlat = TRUE,
  add.demes=TRUE,
  res=600,
  add.grid=TRUE,
  col.demes="black",
  col.grid = "gray20",
  m.plot.xy = {
    map("state", regions = c("california", "oregon"), col = "black", lwd = 2, add = TRUE)
  },
  q.plot.xy = {
    map("state", regions = c("california", "oregon"), col = "black", lwd = 2, add = TRUE)
  }
)















































# Check if the rEEMSplots source directory exists
if [ -d "./rEEMSplots" ]; then
  # Install rEEMSplots package from source
  R CMD INSTALL ./rEEMSplots
else
  echo "Error: Move to the directory that contains the rEEMSplots source to install the package."
  exit 1
fi

# Define the paths for generating graphics
mcmcpath="./data/barrier-schemeZ-nIndiv300-nSites3000-EEMS-nDemes200-simno1"
plotpath="./plot/barrier-schemeZ-nIndiv300-nSites3000-EEMS-nDemes200-simno1-rEEMSplots"

# Run R script to generate graphics
Rscript - <<EOF
# Load rEEMSplots library
library(rEEMSplots)

# Generate graphics
eems.plots("$mcmcpath", "$plotpath", longlat = TRUE)

# Load objects from the RData file
load(paste0("$plotpath", "-rdist.RData"))

# List objects loaded
ls()

# Reproduce plotpath-rdist01.png
library(ggplot2)
library(dplyr)
ggplot(B.component %>% filter(size > 1),
       aes(fitted, obsrvd)) +
  geom_point()
EOF

# Optionally, you might want to do something with the output files, like move them to a different directory
# mv /path/to/output_directory/* /path/to/destination/directory/
