library(yhat)
library(ggplot2)
library(boot)

out_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/09_mmrr"

genMat <- as.matrix(read.table(paste0(out_dir, "/gendist.all.txt"), header = FALSE))
geoMat <- as.matrix(read.table(paste0(out_dir, "/geodist.all.txt"), header = FALSE))
ecoMat <- as.matrix(read.table(paste0(out_dir, "/ecodist.All.txt"), header = FALSE))
pc1Mat <- as.matrix(read.table(paste0(out_dir, "/pc1dist.all.txt"), header = FALSE))

cat("genMat:", dim(genMat), "\n")
cat("geoMat:", dim(geoMat), "\n")
cat("ecoMat:", dim(ecoMat), "\n")
cat("pc1Mat:", dim(pc1Mat), "\n")

if (!all(dim(genMat) == dim(geoMat)) ||
    !all(dim(genMat) == dim(ecoMat)) ||
    !all(dim(genMat) == dim(pc1Mat))) {
  stop("Matrix dimensions do not all match.")
}

genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
pc1Mat_z <- scale(pc1Mat, center = TRUE, scale = TRUE)

unfold <- function(X) {
  x <- vector()
  for (i in 2:nrow(X)) {
    x <- c(x, X[i, 1:(i - 1)])
  }
  return(x)
}

genVec <- unfold(genMat_z)
geoVec <- unfold(geoMat_z)
ecoVec <- unfold(ecoMat_z)
pc1Vec <- unfold(pc1Mat_z)

MMRR <- function(Y, X, nperm = 999) {
  nrowsY <- nrow(Y)
  y <- unfold(Y)

  if (is.null(names(X))) names(X) <- paste("X", 1:length(X), sep = "")

  Xmats <- sapply(X, unfold)
  fit <- lm(y ~ Xmats)
  coeffs <- fit$coefficients
  summ <- summary(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  tprob <- rep(1, length(tstat))
  Fprob <- 1

  for (i in 1:nperm) {
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm)
    fit <- lm(yperm ~ Xmats)
    summ <- summary(fit)
    Fprob <- Fprob + as.numeric(summ$fstatistic[1] >= Fstat)
    tprob <- tprob + as.numeric(abs(summ$coefficients[, "t value"]) >= abs(tstat))
  }

  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  names(r.squared) <- "r.squared"
  names(coeffs) <- c("Intercept", names(X))
  names(tstat) <- paste(c("Intercept", names(X)), "(t)", sep = "")
  names(tp) <- paste(c("Intercept", names(X)), "(p)", sep = "")
  names(Fstat) <- "F-statistic"
  names(Fp) <- "F p-value"

  return(list(
    r.squared = r.squared,
    coefficients = coeffs,
    tstatistic = tstat,
    tpvalue = tp,
    Fstatistic = Fstat,
    Fpvalue = Fp
  ))
}

set.seed(43)

Xmats <- list(
  environment = ecoMat_z,
  structure = pc1Mat_z,
  geography = geoMat_z
)

mmrr_results <- MMRR(genMat_z, Xmats, nperm = 999)
print(mmrr_results)

mmrr_table <- data.frame(
  term = names(mmrr_results$coefficients),
  coefficient = as.numeric(mmrr_results$coefficients),
  t_statistic = as.numeric(mmrr_results$tstatistic),
  p_value = as.numeric(mmrr_results$tpvalue),
  r_squared = as.numeric(mmrr_results$r.squared),
  f_statistic = as.numeric(mmrr_results$Fstatistic),
  f_pvalue = as.numeric(mmrr_results$Fpvalue),
  n_samples = nrow(genMat),
  n_pairwise = length(genVec)
)

write.table(
  mmrr_table,
  paste0(out_dir, "/mmrr_results_all_old_setup.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

data_df <- data.frame(
  GeneticDistance = genVec,
  Geography = geoVec,
  Environment = ecoVec,
  Structure = pc1Vec
)

commonality_results <- commonalityCoefficients(
  data_df,
  "GeneticDistance",
  list("Geography", "Environment", "Structure")
)

print(commonality_results)

commonality_table <- as.data.frame(commonality_results$CC)
commonality_table$Component <- rownames(commonality_table)

write.table(
  commonality_table,
  paste0(out_dir, "/commonality_results_all_old_setup.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

bootstrap_commonality <- function(data, indices) {
  boot_data <- data[indices, ]
  boot_commonality <- commonalityCoefficients(
    boot_data,
    "GeneticDistance",
    list("Geography", "Environment", "Structure")
  )
  return(as.numeric(boot_commonality$CC[, "Coefficient"]))
}

set.seed(123)
boot_results <- boot(data = data_df, statistic = bootstrap_commonality, R = 1000)
ci_results <- apply(boot_results$t, 2, quantile, probs = c(0.025, 0.975))

commonality_ci_table <- data.frame(
  Component = rownames(as.data.frame(commonality_results$CC)),
  Coefficient = as.data.frame(commonality_results$CC)$Coefficient,
  Lower95CI = ci_results[1, ],
  Upper95CI = ci_results[2, ]
)

write.table(
  commonality_ci_table,
  paste0(out_dir, "/commonality_results_all_old_setup_bootstrap_CI.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sink(paste0(out_dir, "/mmrr_commonality_all_old_setup_summary.txt"))
print(mmrr_results)
print(commonality_results)
print(commonality_ci_table)
sink()
