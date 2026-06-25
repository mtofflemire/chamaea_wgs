library(vcfR)
library(adegenet)

# Read VCF data
vcf <- read.vcfR("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/data/ccgp/qc/41-Chamaea.pruned.vcf.gz")
gl1 <- vcfR2genlight(vcf)

# Read metadata
ssw_meta <- read.csv("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/data/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv", header = TRUE)
ssw_meta2 <- ssw_meta[, c("sampleID", "Ecoregion")]
ssw_meta2 <- ssw_meta2[match(gl1$ind.names, ssw_meta2$sampleID), ]


# Convert Ecoregion to factor
ssw_meta2$Ecoregion <- as.factor(ssw_meta2$Ecoregion)



# Assign Ecoregion info to the pop slot of gl1
pop(gl1) <- ssw_meta2$Ecoregion





# Find number of clusters from your dataset
grp <- find.clusters(
  gl1,
  max.n.clust = 3,
  n.pca = 150,
  choose.n.clust = FALSE,
  criterion = "min"
)




find.c
#save BIC plot
pdf(file.path(out_dir, "DAPC_BIC_plot.pdf"), width = 8, height = 6)
plot(seq_along(grp$Kstat), grp$Kstat, type = "b", pch = 19, xlab = "Number of Clusters (K)", ylab = "BIC")
dev.off()








# Which Ecoregions correspond with inferred groups?
pdf(file.path(out_dir, "DAPC_ecoregion_vs_clusters.pdf"), width = 10, height = 8)
table.value(
  table(pop(gl1), grp$grp),
  col.lab = paste("inf", sort(unique(grp$grp))),
  row.lab = levels(pop(gl1))
)
dev.off()




# Run DAPC based on inferred groups
dapc1 <- dapc(gl1, grp$grp, pca.select = "percVar", perc.pca = 80, n.da = length(unique(grp$grp)) - 1)




# Define colors for inferred clusters
grp_levels <- levels(as.factor(grp$grp))
distinct_colors <- rainbow(length(unique(grp$grp)))
dapc_colors <- distinct_colors[seq_along(grp_levels)]




# Basic DAPC scatterplot
pdf(file.path(out_dir, "dapc_plot.pdf"), width = 10, height = 8)
scatter.dapc(dapc1,  cstar = 10, cellipse = 7, xax = 1, ylab = "Discriminant Function 2", legend=T)
dev.off()










