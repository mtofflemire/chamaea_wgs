# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)

MMRR<-function(Y,X,nperm=999){
	#compute regression coefficients and test statistics
	nrowsY<-nrow(Y)
	y<-unfold(Y)
	if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
        Xmats<-sapply(X,unfold)
        fit<-lm(y~Xmats)
	coeffs<-fit$coefficients
	summ<-summary(fit)
	r.squared<-summ$r.squared
	tstat<-summ$coefficients[,"t value"]
	Fstat<-summ$fstatistic[1]
	tprob<-rep(1,length(tstat))
	Fprob<-1

	#perform permutations
	for(i in 1:nperm){
		rand<-sample(1:nrowsY)
		Yperm<-Y[rand,rand]
		yperm<-unfold(Yperm)
		fit<-lm(yperm~Xmats)
		summ<-summary(fit)
                Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
                tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
	}

	#return values
	tp<-tprob/(nperm+1)
	Fp<-Fprob/(nperm+1)
	names(r.squared)<-"r.squared"
	names(coeffs)<-c("Intercept",names(X))
	names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
	names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
	names(Fstat)<-"F-statistic"
	names(Fp)<-"F p-value"
	return(list(r.squared=r.squared,
		coefficients=coeffs,
		tstatistic=tstat,
		tpvalue=tp,
		Fstatistic=Fstat,
		Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
	x<-vector()
	for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
	x<-scale(x, center=TRUE, scale=TRUE)  # Comment this line out if you wish to perform the analysis without standardizing the distance matrices! 
	return(x)
}


# Tutorial for data files gendist.txt, geodist.txt, and ecodist.txt

# Read the matrices from files.
# The read.matrix function requires {tseries} package to be installed and loaded.
# If the files have a row as a header (e.g. column names), then specify 'header=TRUE', default is 'header=FALSE'.
# Load required packages
library(tseries)





# Load distance matrices
genMat <- read.matrix('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/gendist.txt')
geoMat <- read.matrix('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/geodist.txt')
ecoMat <- read.matrix('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/ecodist.txt')
pc1_mat <- read.matrix('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/pc1dist.txt')
  








# Make a list of the explanatory (X) matrices.
# Names are optional.  Order doesn't matter.
# Can include more than two matrices, if desired.
Xmats <- list(
  geographic = geoMat,
  ecology    = ecoMat
)


# Run MMRR function using genMat as the response variable and Xmats as the explanatory variables.
# nperm does not need to be specified, default is nperm=999)
MMRR(genMat,Xmats,nperm=999)

# These data should generate results of approximately:
# Coefficient of geography = 0.778 (p = 0.001)
# Coefficient of ecology = 0.167 (p = 0.063)
# Model r.squared = 0.727 (p = 0.001)
# Note that significance values may change slightly due to the permutation procedure.



y   <- unfold(genMat)
geo <- unfold(geoMat)
eco <- unfold(ecoMat)

res <- MMRR(genMat, Xmats, nperm = 999)

b_geo <- res$coefficients["geographic"]
b_eco <- res$coefficients["ecology"]

y_hat <- b_geo * geo + b_eco * eco

library(ggplot2)

df <- data.frame(
  Observed  = y,
  Predicted = y_hat
)

ggplot(df, aes(Predicted, Observed)) +
  geom_point(alpha = 0.1, color = "gray40") +
  geom_smooth(method = "lm", color = "black", fill = "orange") +
  theme_bw() +
  labs(
    x = "Predicted genetic distance (IBD + IBE)",
    y = "Observed genetic distance"
  )












df_plot <- data.frame(
  GeneticDistance = unfold(genMat),
  EnvironmentalDistance = unfold(ecoMat)
)



ggplot(df_plot, aes(x = EnvironmentalDistance, y = GeneticDistance)) +
  geom_point(alpha = 0.6, color = "gray", size = 1) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Environmental Distance (BIO12)",
    y = "Genetic Distance"
  ) +
  theme_bw()




library(ggplot2)

df_geo <- data.frame(
  GeneticDistance   = unfold(genMat),
  GeographicDistance = unfold(geoMat)
)

ggplot(df_geo, aes(x = GeographicDistance, y = GeneticDistance)) +
  geom_point(alpha = 0.05, size = 0.7, color = "gray") +
  geom_smooth(
    method = "lm",
    color = "black",
    fill = "orange",
    se = TRUE,
    linewidth = 1.2,
    alpha = 0.35
  ) +
  labs(
    x = "Geographic Distance (standardized)",
    y = "Genetic Distance (standardized)"
  ) +
  theme_bw()





df_env <- data.frame(
  GeneticDistance      = unfold(genMat),
  EnvironmentalDistance = unfold(ecoMat)
)

ggplot(df_env, aes(x = EnvironmentalDistance, y = GeneticDistance)) +
  geom_point(alpha = 0.05, size = 0.7, color = "gray") +
  geom_smooth(
    method = "lm",
    color = "black",
    fill = "orange",
    se = TRUE,
    linewidth = 1.2,
    alpha = 0.35
  ) +
  labs(
    x = "Environmental Distance (BIO4: Temperature Seasonality, standardized)",
    y = "Genetic Distance (standardized)"
  ) +
  theme_bw()




# Same pairwise standardized vectors MMRR used
g   <- unfold(genMat)
geo <- unfold(geoMat)
eco <- unfold(ecoMat)

# Remove geographic effect from BOTH variables
g_resid   <- resid(lm(g ~ geo))
eco_resid <- resid(lm(eco ~ geo))

df_partial <- data.frame(
  Environmental_resid = eco_resid,
  Genetic_resid       = g_resid
)

ggplot(df_partial, aes(x = Environmental_resid, y = Genetic_resid)) +
  geom_point(alpha = 0.05, size = 0.7, color = "gray") +
  geom_smooth(
    method = "lm",
    color = "black",
    fill = "orange",
    se = TRUE,
    linewidth = 1.2,
    alpha = 0.35
  ) +
  labs(
    x = "Environmental Distance (BIO4) after removing Geographic Distance",
    y = "Genetic Distance after removing Geographic Distance"
  ) +
  theme_bw()
















library(ggplot2)

# Convert to a long data frame for ggplot2
df_plot <- data.frame(
  GeneticDistance = as.vector(genMat),
  PC1Distance = as.vector(geoMat)
)

ggplot(df_plot, aes(x = ecology, y = GeneticDistance)) +
  geom_point(alpha = 0.1, color = "gray", size = 1) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Geographic Distance (m)",
    y = "Genetic Distance"
  ) +
  theme_bw()

df_plot <- data.frame(
  GeneticDistance = unfold(genMat),
  EnvironmentalDistance = unfold(ecoMat)
)

ggplot(df_plot, aes(x = EnvironmentalDistance, y = GeneticDistance)) +
  geom_point(alpha = 0.1, color = "gray", size = 1) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Environmental Distance (BIO4)",
    y = "Genetic Distance"
  ) +
  theme_bw()



















install.packages("yhat")
library(yhat)
# Load yhat package
library(yhat)

# Convert matrices into vectors
genVec <- unfold(genMat)
geoVec <- unfold(geoMat)
ecoVec <- unfold(ecoMat)
pc1Vec <- unfold(pc1Mat)
unfold <- function(X){
  x <- vector()
  for(i in 2:nrow(X)) x <- c(x, X[i, 1:(i-1)])
  x <- scale(x, center=TRUE, scale=TRUE)  # Standardize the values
  return(x)
}


genVec <- unfold(genMat)
geoVec <- unfold(geoMat)
ecoVec <- unfold(ecoMat)
pc1Vec <- unfold(pc1Mat)



lm_model <- lm(genVec ~ geoVec + ecoVec + pc1Vec)



commonality_results <- commonalityCoefficients(lm_model)

?yhat






?commonalityCoefficients





# Load required package
library(yhat)

# Convert distance matrices into vectors
genVec <- unfold(genMat)
geoVec <- unfold(geoMat)
ecoVec <- unfold(ecoMat)
pc1Vec <- unfold(pc1Mat)

# Create a data frame
data_df <- data.frame(
  GeneticDistance = genVec,
  Geography = geoVec,
  Environment = ecoVec,
  Structure = pc1Vec
)


commonality_results <- commonalityCoefficients(data_df, "GeneticDistance", 
                                               list("Geography", "Environment", "Structure"))


print(commonality_results)


































# Load required packages
install.packages("yhat") # Only run if yhat is not installed
library(yhat)
library(ggplot2)
library(tseries) # Required for read.matrix()

### --- 1. Load and Z-transform Distance Matrices --- ###
# Read matrices from files
genMat <- read.matrix('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/gendist.north.txt')
geoMat <- read.matrix('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/north.geodist.txt')
ecoMat <- read.matrix('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/ecodist.North.txt')
pc1Mat <- read.matrix('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/pc1dist_North.txt')

# Z-transform each matrix before unfolding
genMat_z <- scale(genMat, center=TRUE, scale=TRUE)
geoMat_z <- scale(geoMat, center=TRUE, scale=TRUE)
ecoMat_z <- scale(ecoMat, center=TRUE, scale=TRUE)
pc1Mat_z <- scale(pc1Mat, center=TRUE, scale=TRUE)

# Define an unfold function that converts lower diagonal elements to a vector
unfold <- function(X) {
  x <- vector()
  for (i in 2:nrow(X)) {
    x <- c(x, X[i, 1:(i-1)])
  }
  return(x)  # No need to scale here since we already did it above
}

# Convert matrices to vectors
genVec <- unfold(genMat_z)
geoVec <- unfold(geoMat_z)
ecoVec <- unfold(ecoMat_z)
pc1Vec <- unfold(pc1Mat_z)

### --- 2. Multiple Matrix Regression with Randomization (MMRR) --- ###
MMRR <- function(Y, X, nperm=999) {
  # Compute regression coefficients and test statistics
  nrowsY <- nrow(Y)
  y <- unfold(Y)
  if (is.null(names(X))) names(X) <- paste("X", 1:length(X), sep="")
  Xmats <- sapply(X, unfold)
  fit <- lm(y ~ Xmats)
  coeffs <- fit$coefficients
  summ <- summary(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  tprob <- rep(1, length(tstat))
  Fprob <- 1
  
  # Perform permutations
  for (i in 1:nperm) {
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm)
    fit <- lm(yperm ~ Xmats)
    summ <- summary(fit)
    Fprob <- Fprob + as.numeric(summ$fstatistic[1] >= Fstat)
    tprob <- tprob + as.numeric(abs(summ$coefficients[, "t value"]) >= abs(tstat))
  }
  
  # Return values
  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  names(r.squared) <- "r.squared"
  names(coeffs) <- c("Intercept", names(X))
  names(tstat) <- paste(c("Intercept", names(X)), "(t)", sep="")
  names(tp) <- paste(c("Intercept", names(X)), "(p)", sep="")
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

# Run MMRR
Xmats <- list(environment=ecoMat_z, structure=pc1Mat_z, geography=geoMat_z)
mmrr_results <- MMRR(genMat_z, Xmats, nperm=999)
print(mmrr_results)

### --- 3. Commonality Analysis --- ###
# Create a data frame
data_df <- data.frame(
  GeneticDistance = genVec,
  Geography = geoVec,
  Environment = ecoVec,
  Structure = pc1Vec
)

# Run commonality analysis
commonality_results <- commonalityCoefficients(data_df, "GeneticDistance", 
                                               list("Geography", "Environment", "Structure"))


commona
# Print commonality results
print(commonality_results)

### --- 4. Visualization (Optional) --- ###
df_plot <- data.frame(
  GeneticDistance = genVec,
  PC1Distance = geoVec  # If you meant PC1Mat, swap to pc1Vec
)

ggplot(df_plot, aes(x = PC1Distance, y = GeneticDistance)) +
  geom_point(alpha = 0.3, color = "black", size=0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE, size=0.5) +
  labs(title = "Genetic Distance vs. PC1 Distance",
       x = "PC1 Distance (Genetic Structure)",
       y = "Genetic Distance") +
  theme_minimal()






library(ggplot2)

# Prepare data for plotting
commonality_df <- data.frame(
  Component = c(
    "Unique to Geography", "Unique to Environment", "Unique to Structure",
    "Common Geography & Environment", "Common Geography & Structure", 
    "Common Environment & Structure", "Common Geography, Environment & Structure"
  ),
  Variance = c(0.0003, 0.0219, 0.3906, 0.0086, 0.3019, -0.0151, 0.1583)
)

# Create bar plot
ggplot(commonality_df, aes(x = reorder(Component, Variance), y = Variance, fill = Component)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip for better readability
  theme_minimal() +
  labs(title = "Commonality Analysis Breakdown",
       x = "Variance Component",
       y = "Proportion of Variance Explained") +
  theme(legend.position = "none")












install.packages("boot")  # Install if not already installed
library(boot)
library(yhat)





# Function for bootstrapping commonality coefficients
bootstrap_commonality <- function(data, indices) {
  # Resample data
  boot_data <- data[indices, ]
  
  # Compute commonality coefficients
  commonality_results <- commonalityCoefficients(boot_data, "GeneticDistance", 
                                                 list("Geography", "Environment", "Structure"))
  
  # Extract the coefficients correctly as a numeric vector
  return(as.numeric(commonality_results$CC[, "Coefficient"]))  # Extract as a vector
}





# Prepare data for bootstrapping
data_df <- data.frame(
  GeneticDistance = genVec,
  Geography = geoVec,
  Environment = ecoVec,
  Structure = pc1Vec
)

# Perform bootstrap with 1000 replicates
set.seed(123)  # For reproducibility
boot_results <- boot(data=data_df, statistic=bootstrap_commonality, R=1000)

# Compute 95% confidence intervals using percentile method
ci_results <- boot.ci(boot_results, type="perc")

# Print confidence intervals
print(ci_results)




# Combine CI results with coefficient estimates
# Ensure that CC is a data frame
commonality_cc_df <- as.data.frame(commonality_results$CC)

# Compute confidence intervals using the correct method
ci_results <- apply(boot_results$t, 2, quantile, probs = c(0.025, 0.975))  # 95% CI

# Ensure row names exist for matching components
if (is.null(rownames(commonality_cc_df))) {
  rownames(commonality_cc_df) <- paste0("Component", seq_len(nrow(commonality_cc_df)))
}




# Combine CI results with coefficient estimates
commonality_ci_df <- data.frame(
  Component = rownames(commonality_cc_df),
  Coefficient = commonality_cc_df$Coefficient,  # Correctly extracting coefficients
  Lower95CI = ci_results[1, ],  # 2.5th percentile
  Upper95CI = ci_results[2, ]   # 97.5th percentile
)

# Print results
print(commonality_ci_df)







library(ggplot2)

# Ensure correct ordering (largest to smallest effect)
commonality_ci_df$Component <- factor(commonality_ci_df$Component, 
                                      levels = rev(commonality_ci_df$Component))  # Reverse order for readability

# Create bar plot with error bars
ggplot(commonality_ci_df, aes(x = Coefficient, y = Component, fill = Coefficient > 0)) +
  geom_col(width = 0.5, color = "black") +  # Bars for coefficients
  geom_errorbarh(aes(xmin = Lower95CI, xmax = Upper95CI), height = 0.1, color = "black") +  # Horizontal error bars
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical reference line at 0
  theme_minimal() +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +  # Color bars
  labs(title = "Commonality Coefficients with 95% Confidence Intervals",
       x = "Correlation Coefficient",
       y = NULL) +  # Remove default y-axis label
  theme(
    text = element_text(size = 12),  
    axis.text.y = element_text(size = 10, hjust = 0, vjust = 0.5, margin = margin(r = -10)),  # Push labels left
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    panel.grid.major.y = element_blank(),  # Remove background grid lines
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, -10)  # Remove left margin for tight fit
  ) +
  scale_y_discrete(position = "left")  # Move y-axis labels to the left








# Set up the plot with a standard box and correct axis placement
plot(
  commonality_ci_df$Coefficient,  # X-axis values (coefficients)
  seq_along(commonality_ci_df$Component),  # Y-axis positions
  xlim = range(commonality_ci_df$Lower95CI, commonality_ci_df$Upper95CI),  # Expand limits for CIs
  ylim = c(0.5, nrow(commonality_ci_df) + 0.5),  # Keep spacing readable
  pch = 16,  # Use filled circles for points
  col = ifelse(commonality_ci_df$Coefficient > 0, "blue", "red"),  # Blue for positive, red for negative
  xlab = "Correlation Coefficient",
  ylab = "",
  axes = FALSE  # Remove default axes (we will add custom ones)
)

# Add a standard plot box
box()

# Add error bars (horizontal confidence intervals)
arrows(
  x0 = commonality_ci_df$Lower95CI, 
  x1 = commonality_ci_df$Upper95CI, 
  y0 = seq_along(commonality_ci_df$Component), 
  y1 = seq_along(commonality_ci_df$Component), 
  angle = 90, 
  code = 3, 
  length = 0.05, 
  col = "black"
)

# Add x-axis (bottom) with tick marks
axis(1, tick = TRUE, cex.axis = 0.9)

# Add y-axis on the RIGHT with tick marks and labels
axis(4, at = seq_along(commonality_ci_df$Component), 
     labels = commonality_ci_df$Component, 
     las = 2,  # Rotate text vertically
     tick = TRUE,  # Ensure tick marks appear
     line = 0,  # Align tick marks properly
     cex.axis = 0.9)  # Adjust text size

# Add a vertical reference line at 0
abline(v = 0, lty = 2, col = "black", lwd = 2)

# Add actual percent total values next to each point (right side)
text(
  x = commonality_ci_df$Upper95CI + 0.02,  # Slightly to the right of upper CI
  y = seq_along(commonality_ci_df$Component), 
  labels = paste0(round(commonality_ci_df$`% Total`, 2), "%"),  # Round to 2 decimals
  cex = 0.8,  # Adjust text size
  col = "black",  # Text color
  pos = 4  # Force text to appear to the right
)

# ---- NEW CODE: Add Top Axis for % Total ----
axis(3,  # Adds an x-axis at the top
     at = commonality_ci_df$Coefficient,  # Aligns ticks with coefficient values
     labels = paste0(round(commonality_ci_df$`% Total`, 2), "%"),  # Use % Total as labels
     tick = TRUE,  # Show tick marks
     cex.axis = 0.9,  # Adjust font size
     line = 0.5)  # Moves labels slightly up











library(tidyr)
library(dplyr)

# Name the columns with component labels
colnames(boot_results$t) <- rownames(commonality_results$CC)

# Convert to long format for ggplot
boot_df_long <- as.data.frame(boot_results$t) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "Coefficient")






library(ggplot2)

# Order components by median effect size (optional for readability)
boot_df_long$Component <- factor(
  boot_df_long$Component,
  levels = boot_df_long %>%
    group_by(Component) %>%
    summarize(median = median(Coefficient)) %>%
    arrange(desc(median)) %>%
    pull(Component)
)

# Violin plot
ggplot(boot_df_long, aes(x = Coefficient, y = Component, fill = Component)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Distribution of Bootstrapped Commonality Coefficients",
    x = "Coefficient",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(face = "bold")
  )





ggplot(boot_df_long, aes(x = Coefficient, y = Component, fill = Component)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = "black", linewidth = 0.4) +  # smoother, more readable
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +  # median + IQR inside
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  labs(
    title = "Bootstrapped Distribution of Commonality Coefficients",
    x = "Coefficient",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )








ggplot(boot_df_long, aes(x = Coefficient, y = Component, fill = Component)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  labs(
    title = "Bootstrapped Distribution of Commonality Coefficients",
    x = "Coefficient",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),       # Remove y-axis labels
    axis.ticks.y = element_blank(),      # Remove y-axis ticks
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered title
    plot.margin = margin(10, 10, 10, 10)  # Even margins
  )










library(ggplot2)
library(dplyr)

# Add a new column to position labels to the left of each violin
boot_df_long <- boot_df_long %>%
  mutate(Component = factor(Component, levels = rev(levels(Component))))

label_df <- boot_df_long %>%
  group_by(Component) %>%
  summarize(Center = median(Coefficient))  # Center point for placing label

ggplot(boot_df_long, aes(x = Coefficient, y = Component, fill = Component)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_text(data = label_df, aes(label = Component, x = min(boot_df_long$Coefficient) - 0.01), 
            hjust = 0, color = "black", size = 4) +  # Left-anchored labels
  labs(
    title = "Bootstrapped Distribution of Commonality Coefficients",
    x = "Coefficient",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_cartesian(clip = "off")  # Allow labels outside the plot area









library(ggplot2)
library(dplyr)

# Abbreviations
abbreviations <- c(
  "Pure Geography" = "PG",
  "Pure Environment" = "PE",
  "Pure Structure" = "PS",
  "Geography + Environment" = "G+E",
  "Geography + Structure" = "G+S",
  "Environment + Structure" = "E+S",
  "Geo + Env + Structure" = "G+E+S"
)

# Softer custom color palette
component_colors <- c(
  "PG"     = "#0072B2",  # vibrant blue
  "PE"     = "#009E73",  # vivid green
  "PS"     = "#D55E00",  # strong orange-red
  "G+E"    = "#CC79A7",  # magenta
  "G+S"    = "#F0E442",  # yellow
  "E+S"    = "#56B4E9",  # sky blue
  "G+E+S"  = "#E69F00"   # golden orange
)

# Apply abbreviations
commonality_ci_df$Abbr <- abbreviations[commonality_ci_df$Component]

# Order
commonality_ci_df <- commonality_ci_df %>%
  arrange(desc(Coefficient)) %>%
  mutate(Abbr = factor(Abbr, levels = Abbr))

# Plot
ggplot(commonality_ci_df, aes(x = Coefficient, y = Abbr, fill = Abbr)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbarh(aes(xmin = Lower95CI, xmax = Upper95CI), height = 0.15, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = component_colors, guide = "none") +
  labs(
    title = "Commonality Coefficients with 95% Confidence Intervals",
    x = "Coefficient", y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    panel.grid = element_blank()
  )





component_colors <- c(
  "PG"     = "#1f77b4",
  "PE"     = "#ff7f0e",
  "PS"     = "#2ca02c",
  "G+E"    = "#d62728",
  "G+S"    = "#9467bd",
  "E+S"    = "#8c564b",
  "G+E+S"  = "#17becf"
)


component_colors <- c(
  "PG"     = "#440154",
  "PE"     = "#31688e",
  "PS"     = "#35b779",
  "G+E"    = "#fde725",
  "G+S"    = "#26828e",
  "E+S"    = "#5ec962",
  "G+E+S"  = "#48186a"
)


component_colors <- c(
  "PG"     = "#FF5733",  # fiery red-orange
  "PE"     = "#33C1FF",  # bright sky blue
  "PS"     = "#FF33F6",  # electric magenta
  "G+E"    = "#33FF57",  # vivid green
  "G+S"    = "#FFBD33",  # sunflower yellow-orange
  "E+S"    = "#8D33FF",  # strong violet
  "G+E+S"  = "#33FFF0"   # bright aqua-cyan
)



# Generate 7 blue shades from light to dark
blues <- colorRampPalette(c("#eeeeee", "#555555"))(7)

Gray: c("#eeeeee", "#555555")
# Assign shades to your abbreviations
component_colors <- setNames(blues, c("PG", "PE", "PS", "G+E", "G+S", "E+S", "G+E+S"))
# Add % Total column if missing
# commonality_ci_df$`% Total` <- c(...)  # Your real % values here
commonality_ci_df$Percent <- round(commonality_ci_df$Coefficient / sum(commonality_ci_df$Coefficient) * 100, 2)
ggplot(commonality_ci_df, aes(x = Coefficient, y = Abbr, fill = Abbr)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbarh(aes(xmin = Lower95CI, xmax = Upper95CI), height = 0.15, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_text(
    aes(label = paste0(round(Coefficient * 100, 1), "%"), x = Upper95CI + 0.01),
    hjust = 0, size = 4
  ) +
  scale_fill_manual(values = component_colors, guide = "none") +
  labs(
    x = "Commonality Coefficient",
    y = "Variable Component"  # Y-axis title added here
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )















library(ggplot2)
library(dplyr)

# ---- Example: if your data frame isn't already prepped ----
# Assume commonality_ci_df has these columns:
# Component, Coefficient, Lower95CI, Upper95CI
# Add abbreviations and percent values

abbreviations <- c(
  "Pure Geography" = "PG",
  "Pure Environment" = "PE",
  "Pure Structure" = "PS",
  "Geography + Environment" = "G+E",
  "Geography + Structure" = "G+S",
  "Environment + Structure" = "E+S",
  "Geo + Env + Structure" = "G+E+S"
)

# Add abbrev and order
commonality_ci_df$Abbr <- abbreviations[commonality_ci_df$Component]
commonality_ci_df <- commonality_ci_df %>%
  arrange(desc(Coefficient)) %>%
  mutate(
    Abbr = factor(Abbr, levels = Abbr),
    Percent = round(Coefficient * 100, 1)
  )

# Set all black bars
component_colors <- setNames(rep("black", 7), levels(commonality_ci_df$Abbr))




component_colors <- c(
  "PG"     = "#4B9CD3",  # Carolina blue
  "PE"     = "#E76F51",  # Terracotta
  "PS"     = "#2A9D8F",  # Teal green
  "G+E"    = "#F4A261",  # Apricot
  "G+S"    = "#E9C46A",  # Soft yellow
  "E+S"    = "#A267AC",  # Purple
  "G+E+S"  = "#264653"   # Deep slate
)
# ---- Plot ----
ggplot(commonality_ci_df, aes(x = Coefficient, y = Abbr, fill = Abbr)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbarh(aes(xmin = Lower95CI, xmax = Upper95CI), height = 0.15, color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_text(
    aes(label = paste0(Percent, "%"), x = Upper95CI + 0.01),
    hjust = 0, size = 4
  ) +
  scale_fill_manual(values = component_colors, guide = "none") +
  labs(
    title = "Commonality Coefficients with 95% Confidence Intervals",
    x = "Coefficient",
    y = "Component"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_cartesian(clip = "off")  # So percent labels aren’t clipped

















# Load libraries
library(VennDiagram)
library(grid)
dev.off()
# Open a blank plot
grid.newpage()

# Create the Venn diagram
venn_plot <- draw.triple.venn(
  area1 = 100, area2 = 100, area3 = 100,     # Equal dummy sizes
  n12 = 10, n13 = 10, n23 = 10, n123 = 10,   # Dummy overlaps to force all intersections
  category = c("Geography", "Environment", "Structure"),
  fill = c("yellow", "red", "blue"), # Bright colors
  alpha = 0.5,
  lty = "solid",
  cex = 0,                # Turn off auto labels
  cat.cex = 1.5,
  cat.pos = c(-20, 20, 0),
  cat.dist = 0.08
)
dev.off()
# Manually add your variance partitioning values
# Tweak x/y as needed depending on device size

venn_plot <- draw.triple.venn(
  area1 = 100, area2 = 100, area3 = 100,
  n12 = 10, n13 = 10, n23 = 10, n123 = 10,
  category = c("Geography", "Environment", "Structure"),
  fill = c("blue", "green", "red"),  # blue, green, orange  # <- removes black outlines
  alpha = 0.5,
  lty = "dashed",
  cex = 0,
  cat.cex = 1,
  cat.pos = c(-20, 20, 0),
  cat.dist = 0.08
)
# Pure components
grid.text("0.03%", x = 0.2, y = 0.65)   # Pure Geography
grid.text("2.42%", x = 0.8, y = 0.65)   # Pure Environment
grid.text("45.42%", x = 0.5, y = 0.18)  # Pure Structure

# Pairwise overlaps
grid.text("0.93%", x = 0.5, y = 0.8)    # G ∩ E
grid.text("34.88%", x = 0.32, y = 0.45) # G ∩ S
grid.text("0%", x = 0.68, y = 0.45)     # E ∩ S (was -1.70%)

# 3-way overlap
grid.text("18.02%", x = 0.5, y = 0.6)   # G ∩ E ∩ S

































# Load required packages
library(yhat)
library(ggplot2)
library(tseries) # for read.matrix()

### --- 1. Load and Z-transform Distance Matrices --- ###

genMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/gendist.all.txt"
)

geoMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/geodist.all.txt"
)

ecoMat <- read.matrix(
  '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/ecodist.All.txt'
)

pc1Mat <- read.matrix(
  '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/pc1dist.all.txt'
)

# Check dimensions
cat("genMat:", dim(genMat), "\n")
cat("geoMat:", dim(geoMat), "\n")
cat("ecoMat:", dim(ecoMat), "\n")
cat("pc1Mat:", dim(pc1Mat), "\n")

if (!all(dim(genMat) == dim(geoMat)) ||
    !all(dim(genMat) == dim(ecoMat)) ||
    !all(dim(genMat) == dim(pc1Mat))) {
  stop("Matrix dimensions do not all match.")
}

# Z-transform each matrix before unfolding
genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
pc1Mat_z <- scale(pc1Mat, center = TRUE, scale = TRUE)

# Convert lower triangle of a matrix to a vector
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

### --- 2. Multiple Matrix Regression with Randomization --- ###

MMRR <- function(Y, X, nperm = 999) {
  nrowsY <- nrow(Y)
  y <- unfold(Y)
  
  if (is.null(names(X))) {
    names(X) <- paste("X", 1:length(X), sep = "")
  }
  
  Xmats <- as.data.frame(sapply(X, unfold))
  colnames(Xmats) <- names(X)
  
  fit <- lm(y ~ ., data = Xmats)
  summ <- summary(fit)
  
  coeffs <- coef(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  
  tprob <- rep(1, length(tstat))
  names(tprob) <- names(tstat)
  Fprob <- 1
  
  for (i in 1:nperm) {
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm)
    
    fit_perm <- lm(yperm ~ ., data = Xmats)
    summ_perm <- summary(fit_perm)
    
    perm_tstat <- summ_perm$coefficients[, "t value"]
    
    Fprob <- Fprob + as.numeric(summ_perm$fstatistic[1] >= Fstat)
    
    common_terms <- intersect(names(tstat), names(perm_tstat))
    tprob[common_terms] <- tprob[common_terms] +
      as.numeric(abs(perm_tstat[common_terms]) >= abs(tstat[common_terms]))
  }
  
  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  
  names(r.squared) <- "r.squared"
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

Xmats <- list(
  Environment = ecoMat_z,
  Structure = pc1Mat_z,
  Geography = geoMat_z
)

mmrr_results <- MMRR(genMat_z, Xmats, nperm = 999)
print(mmrr_results)

# Optional diagnostic: shows if one predictor was dropped because of collinearity
Xmats_check <- as.data.frame(sapply(Xmats, unfold))
fit_check <- lm(unfold(genMat_z) ~ ., data = Xmats_check)
print(summary(fit_check))
print(alias(fit_check))

### --- 3. Commonality Analysis --- ###

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

### --- 4. Plot Genetic Distance vs PC1 Distance --- ###

df_plot <- data.frame(
  GeneticDistance = genVec,
  PC1Distance = pc1Vec
)

ggplot(df_plot, aes(x = PC1Distance, y = GeneticDistance)) +
  geom_point(alpha = 0.3, color = "black", size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
  labs(
    title = "Genetic Distance vs. PC1 Distance",
    x = "PC1 Distance (Genetic Structure)",
    y = "Genetic Distance"
  ) +
  theme_minimal()

### --- 5. Plot Commonality Analysis Breakdown --- ###

commonality_df <- as.data.frame(commonality_results$CC)
commonality_df$Component <- rownames(commonality_df)

print(commonality_df)
print(names(commonality_df))

numeric_cols <- names(commonality_df)[sapply(commonality_df, is.numeric)]

coef_col <- numeric_cols[
  grepl("coef|common", numeric_cols, ignore.case = TRUE)
][1]

if (is.na(coef_col)) {
  coef_col <- numeric_cols[1]
}

commonality_plot_df <- data.frame(
  Component = commonality_df$Component,
  Variance = commonality_df[[coef_col]]
)

commonality_plot_df <- commonality_plot_df[
  !grepl("total", commonality_plot_df$Component, ignore.case = TRUE),
]

print(commonality_plot_df)

ggplot(commonality_plot_df, aes(x = reorder(Component, Variance), y = Variance, fill = Component)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Commonality Analysis Breakdown",
    x = "Variance Component",
    y = "Proportion of Variance Explained"
  ) +
  theme(legend.position = "none")




library(stringr)

library(stringr)

commonality_plot_df$Component_wrapped <- str_wrap(
  commonality_plot_df$Component,
  width = 28
)

p_commonality <- ggplot(
  commonality_plot_df,
  aes(x = reorder(Component_wrapped, Variance), y = Variance, fill = Component)
) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Proportion of Variance Explained"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    plot.margin = margin(10, 20, 10, 10)
  )

p_commonality




ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/commonality_breakdown_ALL.pdf",
  plot = p_commonality,
  width = 6,
  height = 7
)






































# Load required packages
library(yhat)
library(ggplot2)
library(tseries) # for read.matrix()
library(stringr)

### --- 1. Load and Z-transform Distance Matrices --- ###

genMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/gendist.north.txt"
)

geoMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/geodist.north.txt"
)

ecoMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/ecodist.North.txt"
)

pc1Mat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/pc1dist.North.txt"
)

# Check dimensions
cat("genMat:", dim(genMat), "\n")
cat("geoMat:", dim(geoMat), "\n")
cat("ecoMat:", dim(ecoMat), "\n")
cat("pc1Mat:", dim(pc1Mat), "\n")

if (!all(dim(genMat) == dim(geoMat)) ||
    !all(dim(genMat) == dim(ecoMat)) ||
    !all(dim(genMat) == dim(pc1Mat))) {
  stop("Matrix dimensions do not all match.")
}

# Z-transform each matrix before unfolding
genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
pc1Mat_z <- scale(pc1Mat, center = TRUE, scale = TRUE)

# Convert lower triangle of a matrix to a vector
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

### --- 2. Multiple Matrix Regression with Randomization --- ###

MMRR <- function(Y, X, nperm = 999) {
  nrowsY <- nrow(Y)
  y <- unfold(Y)
  
  if (is.null(names(X))) {
    names(X) <- paste("X", 1:length(X), sep = "")
  }
  
  Xmats <- as.data.frame(sapply(X, unfold))
  colnames(Xmats) <- names(X)
  
  fit <- lm(y ~ ., data = Xmats)
  summ <- summary(fit)
  
  coeffs <- coef(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  
  tprob <- rep(1, length(tstat))
  names(tprob) <- names(tstat)
  Fprob <- 1
  
  for (i in 1:nperm) {
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm)
    
    fit_perm <- lm(yperm ~ ., data = Xmats)
    summ_perm <- summary(fit_perm)
    
    perm_tstat <- summ_perm$coefficients[, "t value"]
    
    Fprob <- Fprob + as.numeric(summ_perm$fstatistic[1] >= Fstat)
    
    common_terms <- intersect(names(tstat), names(perm_tstat))
    tprob[common_terms] <- tprob[common_terms] +
      as.numeric(abs(perm_tstat[common_terms]) >= abs(tstat[common_terms]))
  }
  
  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  
  names(r.squared) <- "r.squared"
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

Xmats <- list(
  Environment = ecoMat_z,
  Structure = pc1Mat_z,
  Geography = geoMat_z
)

mmrr_results <- MMRR(genMat_z, Xmats, nperm = 999)
print(mmrr_results)

# Optional diagnostic: shows if one predictor was dropped because of collinearity
Xmats_check <- as.data.frame(sapply(Xmats, unfold))
fit_check <- lm(unfold(genMat_z) ~ ., data = Xmats_check)
print(summary(fit_check))
print(alias(fit_check))

### --- 3. Commonality Analysis --- ###

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

### --- 4. Plot Genetic Distance vs PC1 Distance --- ###

df_plot <- data.frame(
  GeneticDistance = genVec,
  PC1Distance = pc1Vec
)

p_pc1 <- ggplot(df_plot, aes(x = PC1Distance, y = GeneticDistance)) +
  geom_point(alpha = 0.3, color = "black", size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
  labs(
    title = "North: Genetic Distance vs. PC1 Distance",
    x = "PC1 Distance (Genetic Structure)",
    y = "Genetic Distance"
  ) +
  theme_minimal()

p_pc1

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/genetic_vs_pc1_North.pdf",
  plot = p_pc1,
  width = 6,
  height = 5
)

### --- 5. Plot Commonality Analysis Breakdown --- ###

commonality_df <- as.data.frame(commonality_results$CC)
commonality_df$Component <- rownames(commonality_df)

print(commonality_df)
print(names(commonality_df))

numeric_cols <- names(commonality_df)[sapply(commonality_df, is.numeric)]

coef_col <- numeric_cols[
  grepl("coef|common", numeric_cols, ignore.case = TRUE)
][1]

if (is.na(coef_col)) {
  coef_col <- numeric_cols[1]
}

commonality_plot_df <- data.frame(
  Component = commonality_df$Component,
  Variance = commonality_df[[coef_col]]
)

commonality_plot_df <- commonality_plot_df[
  !grepl("total", commonality_plot_df$Component, ignore.case = TRUE),
]

print(commonality_plot_df)

commonality_plot_df$Component_wrapped <- str_wrap(
  commonality_plot_df$Component,
  width = 28
)

p_commonality <- ggplot(
  commonality_plot_df,
  aes(x = reorder(Component_wrapped, Variance), y = Variance, fill = Component)
) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "North: Commonality Analysis Breakdown",
    x = NULL,
    y = "Proportion of Variance Explained"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    plot.margin = margin(10, 20, 10, 10)
  )

p_commonality

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/commonality_breakdown_North.pdf",
  plot = p_commonality,
  width = 6,
  height = 7
)

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/commonality_breakdown_North.png",
  plot = p_commonality,
  width = 6,
  height = 7,
  dpi = 300
)






























# Load required packages
library(yhat)
library(ggplot2)
library(tseries) # for read.matrix()
library(stringr)

### --- 1. Load and Z-transform Distance Matrices --- ###

genMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/gendist.south.txt"
)

geoMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/geodist.south.txt"
)

ecoMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/ecodist.South.txt"
)

pc1Mat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/pc1dist.South.txt"
)

# Check dimensions
cat("genMat:", dim(genMat), "\n")
cat("geoMat:", dim(geoMat), "\n")
cat("ecoMat:", dim(ecoMat), "\n")
cat("pc1Mat:", dim(pc1Mat), "\n")

if (!all(dim(genMat) == dim(geoMat)) ||
    !all(dim(genMat) == dim(ecoMat)) ||
    !all(dim(genMat) == dim(pc1Mat))) {
  stop("Matrix dimensions do not all match.")
}

# Z-transform each matrix before unfolding
genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
pc1Mat_z <- scale(pc1Mat, center = TRUE, scale = TRUE)

# Convert lower triangle of a matrix to a vector
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

### --- 2. Multiple Matrix Regression with Randomization --- ###

MMRR <- function(Y, X, nperm = 999) {
  nrowsY <- nrow(Y)
  y <- unfold(Y)
  
  if (is.null(names(X))) {
    names(X) <- paste("X", 1:length(X), sep = "")
  }
  
  Xmats <- as.data.frame(sapply(X, unfold))
  colnames(Xmats) <- names(X)
  
  fit <- lm(y ~ ., data = Xmats)
  summ <- summary(fit)
  
  coeffs <- coef(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  
  tprob <- rep(1, length(tstat))
  names(tprob) <- names(tstat)
  Fprob <- 1
  
  for (i in 1:nperm) {
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm)
    
    fit_perm <- lm(yperm ~ ., data = Xmats)
    summ_perm <- summary(fit_perm)
    
    perm_tstat <- summ_perm$coefficients[, "t value"]
    
    Fprob <- Fprob + as.numeric(summ_perm$fstatistic[1] >= Fstat)
    
    common_terms <- intersect(names(tstat), names(perm_tstat))
    tprob[common_terms] <- tprob[common_terms] +
      as.numeric(abs(perm_tstat[common_terms]) >= abs(tstat[common_terms]))
  }
  
  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  
  names(r.squared) <- "r.squared"
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

Xmats <- list(
  Environment = ecoMat_z,
  Structure = pc1Mat_z,
  Geography = geoMat_z
)

mmrr_results <- MMRR(genMat_z, Xmats, nperm = 999)
print(mmrr_results)

# Optional diagnostic: shows if one predictor was dropped because of collinearity
Xmats_check <- as.data.frame(sapply(Xmats, unfold))
fit_check <- lm(unfold(genMat_z) ~ ., data = Xmats_check)
print(summary(fit_check))
print(alias(fit_check))

### --- 3. Commonality Analysis --- ###

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

### --- 4. Plot Genetic Distance vs PC1 Distance --- ###

df_plot <- data.frame(
  GeneticDistance = genVec,
  PC1Distance = pc1Vec
)

p_pc1 <- ggplot(df_plot, aes(x = PC1Distance, y = GeneticDistance)) +
  geom_point(alpha = 0.3, color = "black", size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
  labs(
    title = "South: Genetic Distance vs. PC1 Distance",
    x = "PC1 Distance (Genetic Structure)",
    y = "Genetic Distance"
  ) +
  theme_minimal()

p_pc1

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/genetic_vs_pc1_South.pdf",
  plot = p_pc1,
  width = 6,
  height = 5
)

### --- 5. Plot Commonality Analysis Breakdown --- ###

commonality_df <- as.data.frame(commonality_results$CC)
commonality_df$Component <- rownames(commonality_df)

print(commonality_df)
print(names(commonality_df))

numeric_cols <- names(commonality_df)[sapply(commonality_df, is.numeric)]

coef_col <- numeric_cols[
  grepl("coef|common", numeric_cols, ignore.case = TRUE)
][1]

if (is.na(coef_col)) {
  coef_col <- numeric_cols[1]
}

commonality_plot_df <- data.frame(
  Component = commonality_df$Component,
  Variance = commonality_df[[coef_col]]
)

commonality_plot_df <- commonality_plot_df[
  !grepl("total", commonality_plot_df$Component, ignore.case = TRUE),
]

print(commonality_plot_df)

commonality_plot_df$Component_wrapped <- str_wrap(
  commonality_plot_df$Component,
  width = 28
)

p_commonality <- ggplot(
  commonality_plot_df,
  aes(x = reorder(Component_wrapped, Variance), y = Variance, fill = Component)
) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Proportion of Variance Explained"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    plot.margin = margin(10, 20, 10, 10)
  )

p_commonality

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/commonality_breakdown_South.pdf",
  plot = p_commonality,
  width = 6,
  height = 7
)

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/commonality_breakdown_South.png",
  plot = p_commonality,
  width = 6,
  height = 7,
  dpi = 300
)




















# Load required packages
library(yhat)
library(ggplot2)
library(tseries) # for read.matrix()
library(stringr)

### --- 1. Load and Z-transform Distance Matrices --- ###

genMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/gendist.all.txt"
)

geoMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/geodist.all.txt"
)

ecoMat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/ecodist.All.txt"
)

pc1Mat <- read.matrix(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts/pc1dist.all.txt"
)

# Check dimensions
cat("genMat:", dim(genMat), "\n")
cat("geoMat:", dim(geoMat), "\n")
cat("ecoMat:", dim(ecoMat), "\n")
cat("pc1Mat:", dim(pc1Mat), "\n")

if (!all(dim(genMat) == dim(geoMat)) ||
    !all(dim(genMat) == dim(ecoMat)) ||
    !all(dim(genMat) == dim(pc1Mat))) {
  stop("Matrix dimensions do not all match.")
}

# Z-transform each matrix before unfolding
genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
pc1Mat_z <- scale(pc1Mat, center = TRUE, scale = TRUE)

# Convert lower triangle of a matrix to a vector
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

### --- 2. Multiple Matrix Regression with Randomization --- ###

MMRR <- function(Y, X, nperm = 999) {
  nrowsY <- nrow(Y)
  y <- unfold(Y)
  
  if (is.null(names(X))) {
    names(X) <- paste("X", 1:length(X), sep = "")
  }
  
  Xmats <- as.data.frame(sapply(X, unfold))
  colnames(Xmats) <- names(X)
  
  fit <- lm(y ~ ., data = Xmats)
  summ <- summary(fit)
  
  coeffs <- coef(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  
  tprob <- rep(1, length(tstat))
  names(tprob) <- names(tstat)
  Fprob <- 1
  
  for (i in 1:nperm) {
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm)
    
    fit_perm <- lm(yperm ~ ., data = Xmats)
    summ_perm <- summary(fit_perm)
    
    perm_tstat <- summ_perm$coefficients[, "t value"]
    
    Fprob <- Fprob + as.numeric(summ_perm$fstatistic[1] >= Fstat)
    
    common_terms <- intersect(names(tstat), names(perm_tstat))
    tprob[common_terms] <- tprob[common_terms] +
      as.numeric(abs(perm_tstat[common_terms]) >= abs(tstat[common_terms]))
  }
  
  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  
  names(r.squared) <- "r.squared"
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

Xmats <- list(
  Environment = ecoMat_z,
  Structure = pc1Mat_z,
  Geography = geoMat_z
)

mmrr_results <- MMRR(genMat_z, Xmats, nperm = 999)
print(mmrr_results)

# Optional diagnostic: shows if one predictor was dropped because of collinearity
Xmats_check <- as.data.frame(sapply(Xmats, unfold))
fit_check <- lm(unfold(genMat_z) ~ ., data = Xmats_check)
print(summary(fit_check))
print(alias(fit_check))

### --- 3. Commonality Analysis --- ###

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

### --- 4. Plot Genetic Distance vs PC1 Distance --- ###

df_plot <- data.frame(
  GeneticDistance = genVec,
  PC1Distance = pc1Vec
)

p_pc1 <- ggplot(df_plot, aes(x = PC1Distance, y = GeneticDistance)) +
  geom_point(alpha = 0.3, color = "black", size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
  labs(
    title = "All Samples: Genetic Distance vs. PC1 Distance",
    x = "PC1 Distance (Genetic Structure)",
    y = "Genetic Distance"
  ) +
  theme_minimal()

p_pc1

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/genetic_vs_pc1_All.pdf",
  plot = p_pc1,
  width = 6,
  height = 5
)

### --- 5. Plot Commonality Analysis Breakdown --- ###

commonality_df <- as.data.frame(commonality_results$CC)
commonality_df$Component <- rownames(commonality_df)

print(commonality_df)
print(names(commonality_df))

numeric_cols <- names(commonality_df)[sapply(commonality_df, is.numeric)]

coef_col <- numeric_cols[
  grepl("coef|common", numeric_cols, ignore.case = TRUE)
][1]

if (is.na(coef_col)) {
  coef_col <- numeric_cols[1]
}

commonality_plot_df <- data.frame(
  Component = commonality_df$Component,
  Variance = commonality_df[[coef_col]]
)

commonality_plot_df <- commonality_plot_df[
  !grepl("total", commonality_plot_df$Component, ignore.case = TRUE),
]

print(commonality_plot_df)

commonality_plot_df$Component_wrapped <- str_wrap(
  commonality_plot_df$Component,
  width = 28
)

p_commonality <- ggplot(
  commonality_plot_df,
  aes(x = reorder(Component_wrapped, Variance), y = Variance, fill = Component)
) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "All Samples: Commonality Analysis Breakdown",
    x = NULL,
    y = "Proportion of Variance Explained"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    plot.margin = margin(10, 20, 10, 10)
  )

p_commonality

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/commonality_breakdown_All.pdf",
  plot = p_commonality,
  width = 6,
  height = 7
)

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/results/Structure/commonality_breakdown_All.png",
  plot = p_commonality,
  width = 6,
  height = 7,
  dpi = 300
)





























# Load required packages
library(yhat)
library(ggplot2)
library(tseries)
library(stringr)

project <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea"
script_dir <- file.path(project, "scripts")
out_dir <- file.path(project, "results/Structure")

### --- Helper functions --- ###

read_dist <- function(path) {
  as.matrix(read.table(path, header = FALSE))
}

unfold <- function(X) {
  x <- vector()
  for (i in 2:nrow(X)) {
    x <- c(x, X[i, 1:(i - 1)])
  }
  return(x)
}

clean_component_name <- function(x) {
  x <- trimws(x)
  x <- gsub(",", "", x)
  x <- gsub("\\s+", " ", x)
  
  has_geo <- grepl("Geography", x)
  has_env <- grepl("Environment", x)
  has_str <- grepl("Structure", x)
  
  if (grepl("Unique", x) && has_geo) {
    return("Unique to Geography")
  }
  
  if (grepl("Unique", x) && has_env) {
    return("Unique to Environment")
  }
  
  if (grepl("Unique", x) && has_str) {
    return("Unique to Structure")
  }
  
  if (has_geo && has_env && has_str) {
    return("Common to Geography, Environment, and Structure")
  }
  
  if (has_geo && has_env) {
    return("Common to Geography and Environment")
  }
  
  if (has_geo && has_str) {
    return("Common to Geography and Structure")
  }
  
  if (has_env && has_str) {
    return("Common to Environment and Structure")
  }
  
  return(NA)
}

run_commonality_only <- function(group_label, gen_file, geo_file, eco_file, pc1_file) {
  genMat <- read_dist(gen_file)
  geoMat <- read_dist(geo_file)
  ecoMat <- read_dist(eco_file)
  pc1Mat <- read_dist(pc1_file)
  
  cat("\nRunning:", group_label, "\n")
  cat("genMat:", dim(genMat), "\n")
  cat("geoMat:", dim(geoMat), "\n")
  cat("ecoMat:", dim(ecoMat), "\n")
  cat("pc1Mat:", dim(pc1Mat), "\n")
  
  if (!all(dim(genMat) == dim(geoMat)) ||
      !all(dim(genMat) == dim(ecoMat)) ||
      !all(dim(genMat) == dim(pc1Mat))) {
    stop("Matrix dimensions do not all match for ", group_label)
  }
  
  genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
  geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
  ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
  pc1Mat_z <- scale(pc1Mat, center = TRUE, scale = TRUE)
  
  genVec <- unfold(genMat_z)
  geoVec <- unfold(geoMat_z)
  ecoVec <- unfold(ecoMat_z)
  pc1Vec <- unfold(pc1Mat_z)
  
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
  
  commonality_df <- as.data.frame(commonality_results$CC)
  commonality_df$Component_raw <- rownames(commonality_df)
  
  commonality_plot_df <- data.frame(
    Group = group_label,
    Component_raw = commonality_df$Component_raw,
    Component = sapply(commonality_df$Component_raw, clean_component_name),
    Variance = commonality_df$Coefficient,
    PercentTotal = commonality_df$`    % Total`
  )
  
  commonality_plot_df <- commonality_plot_df[
    !grepl("total", commonality_plot_df$Component_raw, ignore.case = TRUE),
  ]
  
  if (any(is.na(commonality_plot_df$Component))) {
    print(commonality_plot_df)
    stop("Some component names could not be cleaned for ", group_label)
  }
  
  return(commonality_plot_df)
}

### --- Run commonality for All, North, and South --- ###

commonality_all <- run_commonality_only(
  group_label = "All",
  gen_file = file.path(script_dir, "gendist.all.txt"),
  geo_file = file.path(script_dir, "geodist.all.txt"),
  eco_file = file.path(script_dir, "ecodist.All.txt"),
  pc1_file = file.path(script_dir, "pc1dist.all.txt")
)

commonality_north <- run_commonality_only(
  group_label = "North",
  gen_file = file.path(script_dir, "gendist.north.txt"),
  geo_file = file.path(script_dir, "geodist.north.txt"),
  eco_file = file.path(script_dir, "ecodist.North.txt"),
  pc1_file = file.path(script_dir, "pc1dist.North.txt")
)

commonality_south <- run_commonality_only(
  group_label = "South",
  gen_file = file.path(script_dir, "gendist.south.txt"),
  geo_file = file.path(script_dir, "geodist.south.txt"),
  eco_file = file.path(script_dir, "ecodist.South.txt"),
  pc1_file = file.path(script_dir, "pc1dist.South.txt")
)

combined_commonality <- rbind(
  commonality_all,
  commonality_north,
  commonality_south
)

### --- Set fixed component order --- ###

component_order <- c(
  "Unique to Geography",
  "Unique to Environment",
  "Unique to Structure",
  "Common to Geography and Environment",
  "Common to Geography and Structure",
  "Common to Environment and Structure",
  "Common to Geography, Environment, and Structure"
)

combined_commonality$Component <- factor(
  combined_commonality$Component,
  levels = rev(component_order)
)

combined_commonality$Component_wrapped <- str_wrap(
  as.character(combined_commonality$Component),
  width = 32
)

combined_commonality$Component_wrapped <- factor(
  combined_commonality$Component_wrapped,
  levels = str_wrap(rev(component_order), width = 32)
)

combined_commonality$Group <- factor(
  combined_commonality$Group,
  levels = c("All", "North", "South")
)

print(combined_commonality)

write.table(
  combined_commonality,
  file.path(out_dir, "commonality_combined_All_North_South.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- Combined grouped bar plot --- ###

p_combined_commonality <- ggplot(
  combined_commonality,
  aes(x = Component_wrapped, y = Variance, fill = Group)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "All" = "green3",
      "North" = "blue",
      "South" = "red"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Commonality Analysis Across All, North, and South",
    x = NULL,
    y = "Proportion of Variance Explained",
    fill = "Sample Set"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.position = "top",
    plot.margin = margin(10, 20, 10, 10)
  )

p_combined_commonality

ggsave(
  filename = file.path(out_dir, "commonality_combined_All_North_South.pdf"),
  plot = p_combined_commonality,
  width = 9,
  height = 7
)

ggsave(
  filename = file.path(out_dir, "commonality_combined_All_North_South.png"),
  plot = p_combined_commonality,
  width = 9,
  height = 7,
  dpi = 300
)


























# Load required packages
library(yhat)
library(ggplot2)
library(tseries)
library(stringr)

project <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea"
script_dir <- file.path(project, "scripts")
out_dir <- file.path(project, "results/Structure")

### --- Helper functions --- ###

read_dist <- function(path) {
  as.matrix(read.table(path, header = FALSE))
}

unfold <- function(X) {
  x <- vector()
  for (i in 2:nrow(X)) {
    x <- c(x, X[i, 1:(i - 1)])
  }
  return(x)
}

clean_component_name <- function(x) {
  x <- trimws(x)
  x <- gsub(",", "", x)
  x <- gsub("\\s+", " ", x)
  
  has_geo <- grepl("Geography", x)
  has_env <- grepl("Environment", x)
  has_str <- grepl("Structure", x)
  
  if (grepl("Unique", x) && has_geo) {
    return("Unique Geography")
  }
  
  if (grepl("Unique", x) && has_env) {
    return("Unique Environment")
  }
  
  if (grepl("Unique", x) && has_str) {
    return("Unique Structure")
  }
  
  if (has_geo && has_env && has_str) {
    return("Shared All Three")
  }
  
  if (has_geo && has_env) {
    return("Shared Geography + Environment")
  }
  
  if (has_geo && has_str) {
    return("Shared Geography + Structure")
  }
  
  if (has_env && has_str) {
    return("Shared Environment + Structure")
  }
  
  return(NA)
}

run_commonality_only <- function(group_label, gen_file, geo_file, eco_file, pc1_file) {
  genMat <- read_dist(gen_file)
  geoMat <- read_dist(geo_file)
  ecoMat <- read_dist(eco_file)
  pc1Mat <- read_dist(pc1_file)
  
  cat("\nRunning:", group_label, "\n")
  cat("genMat:", dim(genMat), "\n")
  cat("geoMat:", dim(geoMat), "\n")
  cat("ecoMat:", dim(ecoMat), "\n")
  cat("pc1Mat:", dim(pc1Mat), "\n")
  
  if (!all(dim(genMat) == dim(geoMat)) ||
      !all(dim(genMat) == dim(ecoMat)) ||
      !all(dim(genMat) == dim(pc1Mat))) {
    stop("Matrix dimensions do not all match for ", group_label)
  }
  
  genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
  geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
  ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
  pc1Mat_z <- scale(pc1Mat, center = TRUE, scale = TRUE)
  
  data_df <- data.frame(
    GeneticDistance = unfold(genMat_z),
    Geography = unfold(geoMat_z),
    Environment = unfold(ecoMat_z),
    Structure = unfold(pc1Mat_z)
  )
  
  commonality_results <- commonalityCoefficients(
    data_df,
    "GeneticDistance",
    list("Geography", "Environment", "Structure")
  )
  
  commonality_df <- as.data.frame(commonality_results$CC)
  commonality_df$Component_raw <- rownames(commonality_df)
  
  commonality_plot_df <- data.frame(
    Group = group_label,
    Component_raw = commonality_df$Component_raw,
    Component = sapply(commonality_df$Component_raw, clean_component_name),
    Coefficient = commonality_df$Coefficient,
    PercentTotal = commonality_df$`    % Total`
  )
  
  commonality_plot_df <- commonality_plot_df[
    !grepl("total", commonality_plot_df$Component_raw, ignore.case = TRUE),
  ]
  
  if (any(is.na(commonality_plot_df$Component))) {
    print(commonality_plot_df)
    stop("Some component names could not be cleaned for ", group_label)
  }
  
  return(commonality_plot_df)
}

### --- Run commonality for All, North, and South --- ###

commonality_all <- run_commonality_only(
  group_label = "All",
  gen_file = file.path(script_dir, "gendist.all.txt"),
  geo_file = file.path(script_dir, "geodist.all.txt"),
  eco_file = file.path(script_dir, "ecodist.All.txt"),
  pc1_file = file.path(script_dir, "pc1dist.all.txt")
)

commonality_north <- run_commonality_only(
  group_label = "North",
  gen_file = file.path(script_dir, "gendist.north.txt"),
  geo_file = file.path(script_dir, "geodist.north.txt"),
  eco_file = file.path(script_dir, "ecodist.North.txt"),
  pc1_file = file.path(script_dir, "pc1dist.North.txt")
)

commonality_south <- run_commonality_only(
  group_label = "South",
  gen_file = file.path(script_dir, "gendist.south.txt"),
  geo_file = file.path(script_dir, "geodist.south.txt"),
  eco_file = file.path(script_dir, "ecodist.South.txt"),
  pc1_file = file.path(script_dir, "pc1dist.South.txt")
)

combined_commonality <- rbind(
  commonality_all,
  commonality_north,
  commonality_south
)

### --- Set order and save table --- ###

component_order <- c(
  "Unique Geography",
  "Unique Environment",
  "Unique Structure",
  "Shared Geography + Environment",
  "Shared Geography + Structure",
  "Shared Environment + Structure",
  "Shared All Three"
)

combined_commonality$Group <- factor(
  combined_commonality$Group,
  levels = c("All", "North", "South")
)

combined_commonality$Component <- factor(
  combined_commonality$Component,
  levels = component_order
)

print(combined_commonality)

write.table(
  combined_commonality,
  file.path(out_dir, "commonality_stacked_All_North_South_clean.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- Stacked percent plot --- ###

p_stacked_percent <- ggplot(
  combined_commonality,
  aes(x = Group, y = PercentTotal, fill = Component)
) +
  geom_bar(
    stat = "identity",
    width = 0.75,
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Unique Geography" = "#56B4E9",
      "Unique Environment" = "#009E73",
      "Unique Structure" = "#CC79A7",
      "Shared Geography + Environment" = "#E69F00",
      "Shared Geography + Structure" = "#0072B2",
      "Shared Environment + Structure" = "#D55E00",
      "Shared All Three" = "#F0E442"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Percent of Explained Variance",
    fill = "Component"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 20, 10, 10)
  )

p_stacked_percent

ggsave(
  filename = file.path(out_dir, "commonality_stacked_percent_All_North_South_clean.pdf"),
  plot = p_stacked_percent,
  width = 9,
  height = 6
)

ggsave(
  filename = file.path(out_dir, "commonality_stacked_percent_All_North_South_clean.png"),
  plot = p_stacked_percent,
  width = 9,
  height = 6,
  dpi = 300
)












































# Load required packages
library(yhat)
library(ggplot2)
library(tseries)
library(stringr)

project <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea"
script_dir <- file.path(project, "scripts")
out_dir <- file.path(project, "results/Structure")

### --- Helper functions --- ###

read_dist <- function(path) {
  as.matrix(read.table(path, header = FALSE))
}

unfold <- function(X) {
  x <- vector()
  for (i in 2:nrow(X)) {
    x <- c(x, X[i, 1:(i - 1)])
  }
  return(x)
}

clean_component_name_two_predictors <- function(x) {
  x <- trimws(x)
  x <- gsub(",", "", x)
  x <- gsub("\\s+", " ", x)
  
  has_geo <- grepl("Geography", x)
  has_env <- grepl("Environment", x)
  
  if (grepl("Unique", x) && has_geo) {
    return("Unique Geography")
  }
  
  if (grepl("Unique", x) && has_env) {
    return("Unique Environment")
  }
  
  if (has_geo && has_env) {
    return("Shared Geography + Environment")
  }
  
  return(NA)
}

run_commonality_no_pc1 <- function(group_label, gen_file, geo_file, eco_file) {
  genMat <- read_dist(gen_file)
  geoMat <- read_dist(geo_file)
  ecoMat <- read_dist(eco_file)
  
  cat("\nRunning:", group_label, "\n")
  cat("genMat:", dim(genMat), "\n")
  cat("geoMat:", dim(geoMat), "\n")
  cat("ecoMat:", dim(ecoMat), "\n")
  
  if (!all(dim(genMat) == dim(geoMat)) ||
      !all(dim(genMat) == dim(ecoMat))) {
    stop("Matrix dimensions do not all match for ", group_label)
  }
  
  genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
  geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
  ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
  
  data_df <- data.frame(
    GeneticDistance = unfold(genMat_z),
    Geography = unfold(geoMat_z),
    Environment = unfold(ecoMat_z)
  )
  
  commonality_results <- commonalityCoefficients(
    data_df,
    "GeneticDistance",
    list("Geography", "Environment")
  )
  
  commonality_df <- as.data.frame(commonality_results$CC)
  commonality_df$Component_raw <- rownames(commonality_df)
  
  commonality_plot_df <- data.frame(
    Group = group_label,
    Component_raw = commonality_df$Component_raw,
    Component = sapply(commonality_df$Component_raw, clean_component_name_two_predictors),
    Coefficient = commonality_df$Coefficient,
    PercentTotal = commonality_df$`    % Total`
  )
  
  commonality_plot_df <- commonality_plot_df[
    !grepl("total", commonality_plot_df$Component_raw, ignore.case = TRUE),
  ]
  
  if (any(is.na(commonality_plot_df$Component))) {
    print(commonality_plot_df)
    stop("Some component names could not be cleaned for ", group_label)
  }
  
  return(commonality_plot_df)
}

### --- Run commonality for All, North, and South, without PC1 --- ###

commonality_all <- run_commonality_no_pc1(
  group_label = "All",
  gen_file = file.path(script_dir, "gendist.all.txt"),
  geo_file = file.path(script_dir, "geodist.all.txt"),
  eco_file = file.path(script_dir, "ecodist.All.txt")
)

commonality_north <- run_commonality_no_pc1(
  group_label = "North",
  gen_file = file.path(script_dir, "gendist.north.txt"),
  geo_file = file.path(script_dir, "geodist.north.txt"),
  eco_file = file.path(script_dir, "ecodist.North.txt")
)

commonality_south <- run_commonality_no_pc1(
  group_label = "South",
  gen_file = file.path(script_dir, "gendist.south.txt"),
  geo_file = file.path(script_dir, "geodist.south.txt"),
  eco_file = file.path(script_dir, "ecodist.South.txt")
)

combined_commonality <- rbind(
  commonality_all,
  commonality_north,
  commonality_south
)

component_order <- c(
  "Unique Geography",
  "Unique Environment",
  "Shared Geography + Environment"
)

combined_commonality$Group <- factor(
  combined_commonality$Group,
  levels = c("All", "North", "South")
)

combined_commonality$Component <- factor(
  combined_commonality$Component,
  levels = component_order
)

print(combined_commonality)

write.table(
  combined_commonality,
  file.path(out_dir, "commonality_noPC1_All_North_South.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- Plot 1: Grouped bar plot, raw coefficients --- ###

combined_commonality$Component_wrapped <- str_wrap(
  as.character(combined_commonality$Component),
  width = 28
)

combined_commonality$Component_wrapped <- factor(
  combined_commonality$Component_wrapped,
  levels = str_wrap(rev(component_order), width = 28)
)

p_grouped_no_pc1 <- ggplot(
  combined_commonality,
  aes(x = Component_wrapped, y = Coefficient, fill = Group)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "All" = "green3",
      "North" = "blue",
      "South" = "red"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Commonality Analysis Without PC1",
    x = NULL,
    y = "Proportion of Variance Explained",
    fill = "Sample Set"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.position = "top",
    plot.margin = margin(10, 20, 10, 10)
  )

p_grouped_no_pc1

ggsave(
  filename = file.path(out_dir, "commonality_grouped_noPC1_All_North_South.pdf"),
  plot = p_grouped_no_pc1,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(out_dir, "commonality_grouped_noPC1_All_North_South.png"),
  plot = p_grouped_no_pc1,
  width = 8,
  height = 5,
  dpi = 300
)

### --- Plot 2: Stacked percent plot --- ###

p_stacked_no_pc1 <- ggplot(
  combined_commonality,
  aes(x = Group, y = PercentTotal, fill = Component)
) +
  geom_bar(
    stat = "identity",
    width = 0.75,
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Unique Geography" = "#56B4E9",
      "Unique Environment" = "#009E73",
      "Shared Geography + Environment" = "#E69F00"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Commonality Components Without PC1",
    x = NULL,
    y = "Percent of Explained Variance",
    fill = "Component"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 20, 10, 10)
  )

p_stacked_no_pc1

ggsave(
  filename = file.path(out_dir, "commonality_stacked_noPC1_All_North_South.pdf"),
  plot = p_stacked_no_pc1,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(out_dir, "commonality_stacked_noPC1_All_North_South.png"),
  plot = p_stacked_no_pc1,
  width = 8,
  height = 5,
  dpi = 300
)


















# Load required packages
library(yhat)
library(ggplot2)
library(tseries)
library(stringr)

project <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea"
script_dir <- file.path(project, "scripts")
out_dir <- file.path(project, "results/Structure")

### --- Helper functions --- ###

read_dist <- function(path) {
  as.matrix(read.table(path, header = FALSE))
}

unfold <- function(X) {
  x <- vector()
  for (i in 2:nrow(X)) {
    x <- c(x, X[i, 1:(i - 1)])
  }
  return(x)
}

clean_component_name <- function(x) {
  x <- trimws(x)
  x <- gsub(",", "", x)
  x <- gsub("\\s+", " ", x)
  
  has_geo <- grepl("Geography", x)
  has_env <- grepl("Environment", x)
  has_str <- grepl("Structure", x)
  
  if (grepl("Unique", x) && has_geo) {
    return("Unique Geography")
  }
  
  if (grepl("Unique", x) && has_env) {
    return("Unique Environment")
  }
  
  if (grepl("Unique", x) && has_str) {
    return("Unique Structure")
  }
  
  if (has_geo && has_env && has_str) {
    return("Shared All Three")
  }
  
  if (has_geo && has_env) {
    return("Shared Geography + Environment")
  }
  
  if (has_geo && has_str) {
    return("Shared Geography + Structure")
  }
  
  if (has_env && has_str) {
    return("Shared Environment + Structure")
  }
  
  return(NA)
}

run_commonality_model <- function(group_label, gen_file, geo_file, eco_file, pc1_file = NULL) {
  genMat <- read_dist(gen_file)
  geoMat <- read_dist(geo_file)
  ecoMat <- read_dist(eco_file)
  
  cat("\nRunning:", group_label, "\n")
  cat("genMat:", dim(genMat), "\n")
  cat("geoMat:", dim(geoMat), "\n")
  cat("ecoMat:", dim(ecoMat), "\n")
  
  if (!all(dim(genMat) == dim(geoMat)) ||
      !all(dim(genMat) == dim(ecoMat))) {
    stop("Matrix dimensions do not all match for ", group_label)
  }
  
  genMat_z <- scale(genMat, center = TRUE, scale = TRUE)
  geoMat_z <- scale(geoMat, center = TRUE, scale = TRUE)
  ecoMat_z <- scale(ecoMat, center = TRUE, scale = TRUE)
  
  if (!is.null(pc1_file)) {
    pc1Mat <- read_dist(pc1_file)
    cat("pc1Mat:", dim(pc1Mat), "\n")
    
    if (!all(dim(genMat) == dim(pc1Mat))) {
      stop("PC1 matrix dimensions do not match for ", group_label)
    }
    
    pc1Mat_z <- scale(pc1Mat, center = TRUE, scale = TRUE)
    
    data_df <- data.frame(
      GeneticDistance = unfold(genMat_z),
      Geography = unfold(geoMat_z),
      Environment = unfold(ecoMat_z),
      Structure = unfold(pc1Mat_z)
    )
    
    commonality_results <- commonalityCoefficients(
      data_df,
      "GeneticDistance",
      list("Geography", "Environment", "Structure")
    )
  } else {
    data_df <- data.frame(
      GeneticDistance = unfold(genMat_z),
      Geography = unfold(geoMat_z),
      Environment = unfold(ecoMat_z)
    )
    
    commonality_results <- commonalityCoefficients(
      data_df,
      "GeneticDistance",
      list("Geography", "Environment")
    )
  }
  
  commonality_df <- as.data.frame(commonality_results$CC)
  commonality_df$Component_raw <- rownames(commonality_df)
  
  commonality_plot_df <- data.frame(
    Group = group_label,
    Component_raw = commonality_df$Component_raw,
    Component = sapply(commonality_df$Component_raw, clean_component_name),
    Coefficient = commonality_df$Coefficient,
    PercentTotal = commonality_df$`    % Total`
  )
  
  commonality_plot_df <- commonality_plot_df[
    !grepl("total", commonality_plot_df$Component_raw, ignore.case = TRUE),
  ]
  
  if (any(is.na(commonality_plot_df$Component))) {
    print(commonality_plot_df)
    stop("Some component names could not be cleaned for ", group_label)
  }
  
  print(commonality_plot_df)
  
  return(commonality_plot_df)
}

### --- Run models --- ###
# All samples: Geography + Environment + Structure/PC1
commonality_all <- run_commonality_model(
  group_label = "All",
  gen_file = file.path(script_dir, "gendist.all.txt"),
  geo_file = file.path(script_dir, "geodist.all.txt"),
  eco_file = file.path(script_dir, "ecodist.All.txt"),
  pc1_file = file.path(script_dir, "pc1dist.all.txt")
)

# North: Geography + Environment only
commonality_north <- run_commonality_model(
  group_label = "North",
  gen_file = file.path(script_dir, "gendist.north.txt"),
  geo_file = file.path(script_dir, "geodist.north.txt"),
  eco_file = file.path(script_dir, "ecodist.North.txt"),
  pc1_file = NULL
)

# South: Geography + Environment only
commonality_south <- run_commonality_model(
  group_label = "South",
  gen_file = file.path(script_dir, "gendist.south.txt"),
  geo_file = file.path(script_dir, "geodist.south.txt"),
  eco_file = file.path(script_dir, "ecodist.South.txt"),
  pc1_file = NULL
)

combined_commonality <- rbind(
  commonality_all,
  commonality_north,
  commonality_south
)

### --- Add zero rows for components not included in North/South models --- ###

component_order <- c(
  "Unique Geography",
  "Unique Environment",
  "Unique Structure",
  "Shared Geography + Environment",
  "Shared Geography + Structure",
  "Shared Environment + Structure",
  "Shared All Three"
)

group_order <- c("All", "North", "South")

plot_grid <- expand.grid(
  Group = group_order,
  Component = component_order,
  stringsAsFactors = FALSE
)

combined_complete <- merge(
  plot_grid,
  combined_commonality,
  by = c("Group", "Component"),
  all.x = TRUE
)

combined_complete$Coefficient[is.na(combined_complete$Coefficient)] <- 0
combined_complete$PercentTotal[is.na(combined_complete$PercentTotal)] <- 0
combined_complete$Component_raw[is.na(combined_complete$Component_raw)] <- "Not in model"

combined_complete$Group <- factor(
  combined_complete$Group,
  levels = group_order
)

combined_complete$Component <- factor(
  combined_complete$Component,
  levels = component_order
)

print(combined_complete)

write.table(
  combined_complete,
  file.path(out_dir, "commonality_All_withPC1_NorthSouth_noPC1.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- Plot 1: grouped coefficient plot --- ###

combined_complete$Component_wrapped <- str_wrap(
  as.character(combined_complete$Component),
  width = 28
)

combined_complete$Component_wrapped <- factor(
  combined_complete$Component_wrapped,
  levels = str_wrap(rev(component_order), width = 28)
)

p_grouped_mixed <- ggplot(
  combined_complete,
  aes(x = Component_wrapped, y = Coefficient, fill = Group)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "All" = "green3",
      "North" = "blue",
      "South" = "red"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Commonality Analysis: All with PC1, North/South without PC1",
    x = NULL,
    y = "Proportion of Variance Explained",
    fill = "Sample Set"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.position = "top",
    plot.margin = margin(10, 20, 10, 10)
  )

p_grouped_mixed

ggsave(
  filename = file.path(out_dir, "commonality_grouped_All_withPC1_NorthSouth_noPC1.pdf"),
  plot = p_grouped_mixed,
  width = 9,
  height = 7
)

ggsave(
  filename = file.path(out_dir, "commonality_grouped_All_withPC1_NorthSouth_noPC1.png"),
  plot = p_grouped_mixed,
  width = 9,
  height = 7,
  dpi = 300
)

### --- Plot 2: stacked percent plot --- ###

p_stacked_mixed <- ggplot(
  combined_complete,
  aes(x = Group, y = PercentTotal, fill = Component)
) +
  geom_bar(
    stat = "identity",
    width = 0.75,
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Unique Geography" = "#56B4E9",
      "Unique Environment" = "#009E73",
      "Unique Structure" = "#CC79A7",
      "Shared Geography + Environment" = "#E69F00",
      "Shared Geography + Structure" = "#0072B2",
      "Shared Environment + Structure" = "#D55E00",
      "Shared All Three" = "#F0E442"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Commonality Components: All with PC1, North/South without PC1",
    x = NULL,
    y = "Percent of Explained Variance",
    fill = "Component"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 20, 10, 10)
  )

p_stacked_mixed

ggsave(
  filename = file.path(out_dir, "commonality_stacked_All_withPC1_NorthSouth_noPC1.pdf"),
  plot = p_stacked_mixed,
  width = 9,
  height = 6
)

ggsave(
  filename = file.path(out_dir, "commonality_stacked_All_withPC1_NorthSouth_noPC1.png"),
  plot = p_stacked_mixed,
  width = 9,
  height = 6,
  dpi = 300
)

