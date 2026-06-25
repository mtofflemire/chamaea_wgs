

#process Weir Fst
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)

# Define file path
fst_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/10kb-windows-results/South_v_North.windowed.weir.fst"

# Load Fst data
fst_data <- read.table(fst_file, header=TRUE)

# Remove missing values
fst_data <- na.omit(fst_data)

# Ensure CHROM is a factor with proper ordering
fst_data$CHROM <- factor(fst_data$CHROM, levels = unique(fst_data$CHROM))

# Create a chromosome-specific offset to avoid SNP overlap
chrom_offsets <- fst_data %>%
  group_by(CHROM) %>%
  summarise(chr_offset = max(BIN_START, na.rm=TRUE)) %>%
  mutate(chr_offset = cumsum(lag(chr_offset, default = 0)))  # Cumulative offset

# Merge offsets with the main data
fst_data <- fst_data %>%
  left_join(chrom_offsets, by="CHROM") %>%
  mutate(CUMULATIVE_POS = BIN_START + chr_offset)


#calculate summary statistics
# Calculate the mean of WEIGHTED_FST
mean_fst <- mean(fst_data$WEIGHTED_FST, na.rm=TRUE)

# Print the median Fst value
print(mean_fst)

# Identify Fst outlier thresholds
fst_threshold_5 <- quantile(fst_data$WEIGHTED_FST, 0.95, na.rm=TRUE)
fst_threshold_1 <- quantile(fst_data$WEIGHTED_FST, 0.99, na.rm=TRUE)

# Save outliers
high_fst_5 <- fst_data %>% filter(WEIGHTED_FST >= fst_threshold_5)
high_fst_1 <- fst_data %>% filter(WEIGHTED_FST >= fst_threshold_1)

write.table(high_fst_5, "high_fst_5percent.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(high_fst_1, "high_fst_1percent.txt", row.names=FALSE, quote=FALSE, sep="\t")

# Manhattan Plot with evenly spaced chromosomes and y-axis starting at 0
ggplot(fst_data, aes(x=CUMULATIVE_POS, y=WEIGHTED_FST, color=as.factor(CHROM))) +
  geom_point(alpha=0.5, size=0.5) +
  geom_hline(yintercept=fst_threshold_5, linetype="dashed", color="red") +
  geom_hline(yintercept=fst_threshold_1, linetype="dashed", color="black") +
  scale_x_continuous(breaks = tapply(fst_data$CUMULATIVE_POS, fst_data$CHROM, median),
                     labels = levels(fst_data$CHROM)) +
  scale_color_manual(values=rep(c("black", "gray"), length.out=length(unique(fst_data$CHROM)))) +
  labs(
       x="Chromosome",
       y="Fst") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90, hjust=1)) +
  ylim(0, max(fst_data$WEIGHTED_FST, na.rm=TRUE))

#######################End code














#Process pi for South population
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)

# Define file path
pi_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/windowed_stats/pi/10kbp/south/run_2025-07-02_22-02-46/chamaea.south.samps.windowed.pi"

# Load the pi data
pi_data <- read.table(pi_file, header=TRUE)

# Remove missing values
pi_data <- na.omit(pi_data)

# Ensure CHROM is a factor with proper ordering
pi_data$CHROM <- factor(pi_data$CHROM, levels = unique(pi_data$CHROM))

# Create a chromosome-specific offset to avoid SNP overlap
chrom_offsets <- pi_data %>%
  group_by(CHROM) %>%
  summarise(chr_offset = max(BIN_START, na.rm=TRUE)) %>%
  mutate(chr_offset = cumsum(lag(chr_offset, default = 0)))  # Cumulative offset

# Merge offsets with the main data
pi_data <- pi_data %>%
  left_join(chrom_offsets, by="CHROM") %>%
  mutate(CUMULATIVE_POS = BIN_START + chr_offset)

# Identify π outlier thresholds
pi_threshold_5 <- quantile(pi_data$PI, 0.95, na.rm=TRUE)
pi_threshold_1 <- quantile(pi_data$PI, 0.99, na.rm=TRUE)

# Save high π regions
high_pi_5 <- pi_data %>% filter(PI >= pi_threshold_5)
high_pi_1 <- pi_data %>% filter(PI >= pi_threshold_1)

#write.table(high_pi_5, "high_pi_5percent.txt", row.names=FALSE, quote=FALSE, sep="\t")
#write.table(high_pi_1, "high_pi_1percent.txt", row.names=FALSE, quote=FALSE, sep="\t")

# Manhattan Plot for π values
ggplot(pi_data, aes(x=CUMULATIVE_POS, y=PI, color=as.factor(CHROM))) +
  geom_point(alpha=0.5, size=0.5) +
  geom_hline(yintercept=pi_threshold_5, linetype="dashed", color="red") +
  geom_hline(yintercept=pi_threshold_1, linetype="dashed", color="black") +
  scale_x_continuous(breaks = tapply(pi_data$CUMULATIVE_POS, pi_data$CHROM, median),
                     labels = levels(pi_data$CHROM)) +
  scale_color_manual(values=rep(c("black", "gray"), length.out=length(unique(pi_data$CHROM)))) +
  labs(
       x="Chromosome",
       y="π South") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90, hjust=1))

# Calculate overall median π value
median_pi <- median(pi_data$PI, na.rm=TRUE)
print(median_pi)

# Calculate median π per chromosome
median_pi_by_chrom <- pi_data %>%
  group_by(CHROM) %>%
  summarise(median_pi = median(PI, na.rm=TRUE))

print(median_pi_by_chrom)
##########End code










# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)

# ========================
# ==== USER INPUT HERE ===
# ========================
# Define file path to any population's π file
pi_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/windowed_stats/pi/10kbp/south/run_2025-07-02_22-02-46/chamaea.south.samps.windowed.pi"

# Extract output directory and population name
output_dir <- dirname(pi_file)
pop_label <- str_extract(basename(pi_file), "(north|south|central|east|west)")  # extend if needed

# ========================
# ==== PROCESSING π ======
# ========================

# Load the π data
pi_data <- read.table(pi_file, header=TRUE) %>% na.omit()

# Format chromosome factor
pi_data$CHROM <- factor(pi_data$CHROM, levels = unique(pi_data$CHROM))

# Create cumulative genomic position
chrom_offsets <- pi_data %>%
  group_by(CHROM) %>%
  summarise(chr_offset = max(BIN_START, na.rm=TRUE)) %>%
  mutate(chr_offset = cumsum(lag(chr_offset, default = 0)))

pi_data <- pi_data %>%
  left_join(chrom_offsets, by="CHROM") %>%
  mutate(CUMULATIVE_POS = BIN_START + chr_offset)

# Calculate outlier thresholds
threshold_5 <- quantile(pi_data$PI, 0.95, na.rm=TRUE)
threshold_1 <- quantile(pi_data$PI, 0.99, na.rm=TRUE)

# Save high π regions
write.table(
  pi_data %>% filter(PI >= threshold_5),
  file.path(output_dir, paste0("high_pi_", pop_label, "_5percent.txt")),
  row.names=FALSE, quote=FALSE, sep="\t"
)
write.table(
  pi_data %>% filter(PI >= threshold_1),
  file.path(output_dir, paste0("high_pi_", pop_label, "_1percent.txt")),
  row.names=FALSE, quote=FALSE, sep="\t"
)

# ========================
# ==== PLOT SAVING =======
# ========================
png(
  filename = file.path(output_dir, paste0("pi_", pop_label, "_manhattan.png")),
  width = 3000, height = 1800, res = 300
)
ggplot(pi_data, aes(x=CUMULATIVE_POS, y=PI, color=as.factor(CHROM))) +
  geom_point(alpha=0.5, size=0.5) +
  geom_hline(yintercept=threshold_5, linetype="dashed", color="red") +
  geom_hline(yintercept=threshold_1, linetype="dashed", color="black") +
  scale_x_continuous(
    breaks = tapply(pi_data$CUMULATIVE_POS, pi_data$CHROM, median),
    labels = levels(pi_data$CHROM)) +
  scale_color_manual(values=rep(c("black", "gray"), length.out=length(unique(pi_data$CHROM)))) +
  labs(x="Chromosome", y=paste("π -", str_to_title(pop_label))) +
  theme_minimal() +
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1))
dev.off()

# ========================
# ==== SUMMARY STATS =====
# ========================
median_pi <- median(pi_data$PI, na.rm=TRUE)
mean_pi <- mean(pi_data$PI, na.rm=TRUE)
median_pi_by_chrom <- pi_data %>%
  group_by(CHROM) %>%
  summarise(
    median_pi = median(PI, na.rm=TRUE),
    mean_pi = mean(PI, na.rm=TRUE)
  )

summary_file <- file.path(output_dir, paste0("pi_", pop_label, "_summary.txt"))
sink(summary_file)
cat("π Summary Statistics -", str_to_title(pop_label), "Population\n")
cat("--------------------------------------------------\n")
cat(paste("Overall median π:", round(median_pi, 6)), "\n")
cat(paste("Overall mean π:", round(mean_pi, 6)), "\n\n")
cat("Median and Mean π per chromosome:\n")
print(median_pi_by_chrom)
sink()







#calculate Tajima's D
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define file paths
tajima_south_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/10kb-windows-results/south_v_north_tajimasD_south.Tajima.D"
tajima_north_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/10kb-windows-results/south_v_north_tajimasD_north.Tajima.D"

# Read Tajima’s D data
tajima_south_data <- read.table(tajima_south_file, header=TRUE)
tajima_north_data <- read.table(tajima_north_file, header=TRUE)

# Set column names if needed
colnames(tajima_south_data) <- c("CHROM", "BIN_START", "N_SNPS", "TajimaD")
colnames(tajima_north_data) <- c("CHROM", "BIN_START", "N_SNPS", "TajimaD")

# Remove missing values
tajima_south_data <- na.omit(tajima_south_data)
tajima_north_data <- na.omit(tajima_north_data)







# Ensure CHROM is a factor with proper ordering
tajima_south_data$CHROM <- factor(tajima_south_data$CHROM, levels=unique(tajima_south_data$CHROM))
tajima_north_data$CHROM <- factor(tajima_north_data$CHROM, levels=unique(tajima_north_data$CHROM))

# Compute offsets to avoid overlap in Manhattan plot
chrom_offsets_south <- tajima_south_data %>%
  group_by(CHROM) %>%
  summarise(chr_offset = max(BIN_START, na.rm=TRUE)) %>%
  mutate(chr_offset = cumsum(lag(chr_offset, default=0)))

chrom_offsets_north <- tajima_north_data %>%
  group_by(CHROM) %>%
  summarise(chr_offset = max(BIN_START, na.rm=TRUE)) %>%
  mutate(chr_offset = cumsum(lag(chr_offset, default=0)))

# Merge offsets with the main data
tajima_south_data <- tajima_south_data %>%
  left_join(chrom_offsets_south, by="CHROM") %>%
  mutate(CUMULATIVE_POS = BIN_START + chr_offset)

tajima_north_data <- tajima_north_data %>%
  left_join(chrom_offsets_north, by="CHROM") %>%
  mutate(CUMULATIVE_POS = BIN_START + chr_offset)







# Identify outlier thresholds for Tajima’s D
tajima_south_threshold_5 <- quantile(tajima_south_data$TajimaD, 0.95, na.rm=TRUE)
tajima_south_threshold_1 <- quantile(tajima_south_data$TajimaD, 0.99, na.rm=TRUE)

tajima_north_threshold_5 <- quantile(tajima_north_data$TajimaD, 0.95, na.rm=TRUE)
tajima_north_threshold_1 <- quantile(tajima_north_data$TajimaD, 0.99, na.rm=TRUE)

# Save high Tajima’s D regions
high_tajima_south_5 <- tajima_south_data %>% filter(TajimaD >= tajima_south_threshold_5)
high_tajima_south_1 <- tajima_south_data %>% filter(TajimaD >= tajima_south_threshold_1)

high_tajima_north_5 <- tajima_north_data %>% filter(TajimaD >= tajima_north_threshold_5)
high_tajima_north_1 <- tajima_north_data %>% filter(TajimaD >= tajima_north_threshold_1)

# Save outlier data to files (optional)
# write.table(high_tajima_south_5, "high_tajima_south_5percent.txt", row.names=FALSE, quote=FALSE, sep="\t")
# write.table(high_tajima_south_1, "high_tajima_south_1percent.txt", row.names=FALSE, quote=FALSE, sep="\t")
# write.table(high_tajima_north_5, "high_tajima_north_5percent.txt", row.names=FALSE, quote=FALSE, sep="\t")
# write.table(high_tajima_north_1, "high_tajima_north_1percent.txt", row.names=FALSE, quote=FALSE, sep="\t")




# Manhattan Plot for Tajima's D (South Population)
ggplot(tajima_south_data, aes(x=CUMULATIVE_POS, y=TajimaD, color=as.factor(CHROM))) +
  geom_point(alpha=0.5, size=0.5) +
  geom_hline(yintercept=tajima_south_threshold_5, linetype="dashed", color="red") +
  geom_hline(yintercept=tajima_south_threshold_1, linetype="dashed", color="black") +
  scale_x_continuous(breaks = tapply(tajima_south_data$CUMULATIVE_POS, tajima_south_data$CHROM, median),
                     labels = levels(tajima_south_data$CHROM)) +
  scale_color_manual(values=rep(c("black", "gray"), length.out=length(unique(tajima_south_data$CHROM)))) +
  labs(
    x="Chromosome",
    y="Tajima’s D - South") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90, hjust=1))

# Manhattan Plot for Tajima's D (North Population)
ggplot(tajima_north_data, aes(x=CUMULATIVE_POS, y=TajimaD, color=as.factor(CHROM))) +
  geom_point(alpha=0.5, size=0.5) +
  geom_hline(yintercept=tajima_north_threshold_5, linetype="dashed", color="red") +
  geom_hline(yintercept=tajima_north_threshold_1, linetype="dashed", color="black") +
  scale_x_continuous(breaks = tapply(tajima_north_data$CUMULATIVE_POS, tajima_north_data$CHROM, median),
                     labels = levels(tajima_north_data$CHROM)) +
  scale_color_manual(values=rep(c("black", "gray"), length.out=length(unique(tajima_north_data$CHROM)))) +
  labs(
    x="Chromosome",
    y="Tajima’s D - North") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90, hjust=1))






# Compute overall median Tajima’s D
median_tajima_south <- median(tajima_south_data$TajimaD, na.rm=TRUE)
median_tajima_north <- median(tajima_north_data$TajimaD, na.rm=TRUE)

print(median_tajima_south)
print(median_tajima_north)

# Compute median Tajima’s D per chromosome
median_tajima_south_by_chrom <- tajima_south_data %>%
  group_by(CHROM) %>%
  summarise(median_tajimaD = median(TajimaD, na.rm=TRUE))

median_tajima_north_by_chrom <- tajima_north_data %>%
  group_by(CHROM) %>%
  summarise(median_tajimaD = median(TajimaD, na.rm=TRUE))

print(median_tajima_south_by_chrom)
print(median_tajima_north_by_chrom)







#Filter for Oregon and just chromosome 6 for Tajima's D
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define file path for Oregon population
tajima_oregon_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/chamaea-10kb-windowed/Oregon_tajimasD.Tajima.D"

# Read Tajima’s D data
tajima_oregon_data <- read.table(tajima_oregon_file, header=TRUE)

# Set column names if needed
colnames(tajima_oregon_data) <- c("CHROM", "BIN_START", "N_SNPS", "TajimaD")

# Remove missing values
tajima_oregon_data <- na.omit(tajima_oregon_data)

# **Keep only chromosome 6 (JARCOQ010000006.1)**
tajima_chr6 <- tajima_oregon_data %>% filter(CHROM == "JARCOQ010000006.1")

# **Determine the last quarter of chromosome 6**
max_pos <- max(tajima_chr6$BIN_START, na.rm=TRUE)
cutoff <- max_pos - (max_pos / 4)  # Last 25% of the chromosome

# **Filter for last quarter**
tajima_chr6_last_quarter <- tajima_chr6 %>% filter(BIN_START >= cutoff)

# Identify outlier thresholds for the last quarter of chromosome 6
tajima_chr6_threshold_5 <- quantile(tajima_chr6_last_quarter$TajimaD, 0.95, na.rm=TRUE)
tajima_chr6_threshold_1 <- quantile(tajima_chr6_last_quarter$TajimaD, 0.99, na.rm=TRUE)

# **Plot Tajima's D for the last quarter of chromosome 6**
ggplot(tajima_chr6_last_quarter, aes(x=BIN_START, y=TajimaD)) +
  geom_line(alpha=0.5, size=0.5, color="blue") +
  geom_hline(yintercept=tajima_chr6_threshold_5, linetype="dashed", color="red") +
  geom_hline(yintercept=tajima_chr6_threshold_1, linetype="dashed", color="black") +
  labs(
    x="Position on Chromosome 6 (Last Quarter)",
    y="Tajima’s D - Oregon (Chromosome 6 - Last Quarter)") +
  theme_minimal()

# Print median Tajima’s D for the last quarter of chromosome 6
median_tajima_chr6_last_quarter <- median(tajima_chr6_last_quarter$TajimaD, na.rm=TRUE)
print(median_tajima_chr6_last_quarter)


#End Code





#plot tajima's D for specific part of genome
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define file path for Oregon population
tajima_oregon_file <- "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/chamaea-10kb-windowed/Oregon_tajimasD.Tajima.D"

# Read Tajima’s D data
tajima_oregon_data <- read.table(tajima_oregon_file, header=TRUE)

# Set column names if needed
colnames(tajima_oregon_data) <- c("CHROM", "BIN_START", "N_SNPS", "TajimaD")

# Remove missing values
tajima_oregon_data <- na.omit(tajima_oregon_data)

# **Keep only chromosome 6 (JARCOQ010000006.1)**
tajima_chr6 <- tajima_oregon_data %>% filter(CHROM == "JARCOQ010000006.1")

# **Define SNP position and zoom-in range**
snp_position <- 69000001
zoom_range_start <- snp_position - 1000000  # 50kb before SNP
zoom_range_end <- snp_position + 1000000  # 50kb after SNP

# **Filter data to only keep positions within this small range**
tajima_chr6_zoomed <- tajima_chr6 %>% filter(BIN_START >= zoom_range_start & BIN_START <= zoom_range_end)

# Identify outlier thresholds for the zoomed region
tajima_chr6_threshold_5 <- quantile(tajima_chr6_zoomed$TajimaD, 0.95, na.rm=TRUE)
tajima_chr6_threshold_1 <- quantile(tajima_chr6_zoomed$TajimaD, 0.99, na.rm=TRUE)

# **Plot Tajima's D for the zoomed-in region**
ggplot(tajima_chr6_zoomed, aes(x=BIN_START, y=TajimaD)) +
  geom_line(alpha=0.5, size=1, color="blue") +  # Larger points for zoomed-in view
  geom_hline(yintercept=tajima_chr6_threshold_5, linetype="dashed", color="red") +
  geom_hline(yintercept=tajima_chr6_threshold_1, linetype="dashed", color="black") +
  geom_vline(xintercept=snp_position, linetype="solid", color="red", linewidth=1) +  # Vertical line at SNP
  annotate("text", x=snp_position, y=max(tajima_chr6_zoomed$TajimaD, na.rm=TRUE), 
           label="SNP Position", color="black", angle=90, vjust=-0.5) +  # Label for SNP
  labs(
    x="Position on Chromosome 6 (Zoomed: ±50kb)",
    y="Tajima’s D - Oregon (Zoomed View)") +
  theme_bw()




# **Plot Tajima's D for the zoomed-in region**
ggplot(tajima_chr6_zoomed, aes(x=BIN_START, y=TajimaD)) +
  geom_line(alpha=1, size=1, color="navyblue") +  # Line plot instead of points
  geom_hline(yintercept=tajima_chr6_threshold_5, linetype="dashed", color="red") +
  geom_hline(yintercept=tajima_chr6_threshold_1, linetype="dashed", color="black") +
  geom_vline(xintercept=snp_position, linetype="solid", color="red", linewidth=1) +  # Vertical line at SNP
  annotate("text", x=snp_position, y=max(tajima_chr6_zoomed$TajimaD, na.rm=TRUE), 
           label="SNP Position", color="black", angle=90, vjust=-0.5) +  # Label for SNP
  labs(
    x="Position on Chromosome 6 (Zoomed: ±50kb)",
    y="Tajima’s D - Oregon (Zoomed View)") +
  theme_bw() +  # Base theme with white background
  theme(
    panel.grid = element_blank(),  # Removes gridlines
    panel.border = element_rect(color="black", size=1.5)  # Makes the box thicker
  )
















# Load necessary libraries
library(dplyr)

# Create a summary table
summary_table <- data.frame(
  Population = c("South", "North"),
  
  # Fst values (Only one overall Fst value)
  Mean_Fst = c(mean_fst, mean_fst),  
  
  # Pi (Nucleotide Diversity)
  Median_Pi = c(median_pi, median_pi_north),
  
  # Tajima's D
  Median_TajimaD = c(median_tajima_south, median_tajima_north)
)

# Print summary table
print(summary_table)


write.csv(summary_table, "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/windowed_analysis_summary.csv", row.names=FALSE)

