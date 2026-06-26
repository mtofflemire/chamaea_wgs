library(vcfR)
library(adegenet)



#define paths
args_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(args_file)) dirname(args_file) else getwd()

vcf_file <- file.path(script_dir, "Chamaea_auto_filteredQC_maf0.05_hwe0.01_prune1kb.vcf.gz")
metadata_file <- file.path(script_dir, "CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv")
out_dir <- script_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



#read vcf
vcf <- read.vcfR(vcf_file)
gl1 <- vcfR2genlight(vcf)



#read metadata
ssw_meta <- read.csv(metadata_file, header = TRUE)
ssw_meta2 <- ssw_meta[, c("sampleID", "Ecoregion")]
ssw_meta2 <- ssw_meta2[match(gl1$ind.names, ssw_meta2$sampleID), ]
ssw_meta2$Ecoregion <- as.factor(ssw_meta2$Ecoregion)

pop(gl1) <- ssw_meta2$Ecoregion



#find clusters
grp <- find.clusters(
  gl1,
  max.n.clust = 10,
  n.pca = 150,
  choose.n.clust = FALSE,
  criterion = "min"
)



#save BIC plot
pdf(file.path(out_dir, "DAPC_BIC_plot.pdf"), width = 8, height = 6)
plot(
  seq_along(grp$Kstat),
  grp$Kstat,
  type = "b",
  pch = 19,
  xlab = "Number of Clusters (K)",
  ylab = "BIC"
)
dev.off()



#which ecoregions correspond with inferred groups?
pdf(file.path(out_dir, "DAPC_ecoregion_vs_clusters.pdf"), width = 10, height = 8)
table.value(
  table(pop(gl1), grp$grp),
  col.lab = paste("inf", sort(unique(grp$grp))),
  row.lab = levels(pop(gl1))
)
dev.off()

write.csv(
  table(pop(gl1), grp$grp),
  file.path(out_dir, "DAPC_ecoregion_vs_clusters.csv")
)



#run dapc
dapc1 <- dapc(
  gl1,
  grp$grp,
  pca.select = "percVar",
  perc.pca = 80,
  n.da = length(unique(grp$grp)) - 1
)

saveRDS(grp, file.path(out_dir, "DAPC_find_clusters.rds"))
saveRDS(dapc1, file.path(out_dir, "DAPC_results.rds"))



#plot dapc
pdf(file.path(out_dir, "DAPC_plot.pdf"), width = 10, height = 8)
scatter.dapc(
  dapc1,
  cstar = 10,
  cellipse = 7,
  xax = 1,
  ylab = "Discriminant Function 2",
  legend = TRUE
)
dev.off()
