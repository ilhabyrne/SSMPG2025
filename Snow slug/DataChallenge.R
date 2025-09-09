# Byrne, 2025
# SSMPG - Data Challenge

# Snow slug simulated dataset

setwd("~/Documents/GitHub/SSMPG2025/Snow slug")


# Population structure ----

# Read your specific data format (semicolon-separated)
geno_data <- read.csv("genotype_frequencies.csv", 
                      sep = ";", 
                      header = TRUE, 
                      row.names = 1)  # Use first column as row names

# Quick data summary
head(geno_data[, 1:10])  # Show first 10 columns
summary(geno_data[, 1:5])  # Summary of first 5 columns

# Handle any missing values if present
geno_data_clean <- na.omit(geno_data)

# Optional: Remove variants with low variation (MAF < 0.05)
variant_means <- colMeans(geno_data_clean, na.rm = TRUE)
variant_vars <- apply(geno_data_clean, 2, var, na.rm = TRUE)

# Keep variants with reasonable variation
high_var_variants <- variant_vars > 0.01  # Adjust threshold as needed
geno_filtered <- geno_data_clean[, high_var_variants]

# Perform PCA
pca_result <- prcomp(geno_data, 
                     center = TRUE, 
                     scale. = TRUE)

# View PCA summary
summary(pca_result)

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2,] * 100
cumvar_explained <- summary(pca_result)$importance[3,] * 100

# Extract PC scores with sample names
pc_scores <- data.frame(pca_result$x)
pc_scores$Sample <- rownames(pc_scores)
pc_scores$Sample_ID <- paste("Sample", 1:nrow(pc_scores))

# Create visualizations
# Scree plot
scree_data <- data.frame(
  PC = 1:min(20, ncol(pca_result$x)),
  Variance = var_explained[1:min(20, ncol(pca_result$x))]
)

p1 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  labs(x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_classic() +
  scale_x_continuous(breaks = 1:20)

print(p1)

# PC1 vs PC2 plot
p2 <- ggplot(pc_scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 4, alpha = 0.8, color = "darkblue") +
  labs(x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")) +
  theme_classic()

print(p2)

# PC1 vs PC3 plot
p3 <- ggplot(pc_scores, aes(x = PC1, y = PC3)) +
  geom_point(size = 4, alpha = 0.8, color = "darkblue") +
  labs(x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
       y = paste0("PC3 (", round(var_explained[3], 1), "% variance)")) +
  theme_classic()

print(p3)

# PC2 vs PC3 plot
p4 <- ggplot(pc_scores, aes(x = PC2, y = PC3)) +
  geom_point(size = 4, alpha = 0.8, color = "darkblue") +
  labs(x = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
       y = paste0("PC3 (", round(var_explained[3], 1), "% variance)")) +
  theme_classic()

print(p4)


# GEA analysis using LFMM2 ----

# Convert to matrix format (samples x loci)
Y <- as.matrix(geno_data)
cat("Y matrix dimensions:", dim(Y), "\n")

# Read your environmental data
env_data <- read.csv("observed_environment.csv",
                     sep = ";",
                     header = TRUE, 
                     row.names = 1)  # Assuming first column has sample names

# Standardize environmental variables
X <- as.matrix(env_data)
X_scaled <- scale(X)
X <- X_scaled

# Check for missing data and handle if necessary
missing_geno <- sum(is.na(Y))
missing_env <- sum(is.na(X))

# Run LFMM2 analysis for each environmental variable
# Initialize results storage
results_list <- list()
pvalue_matrices <- list()

# Determine number of latent factors (K)
# You can use the K from your SNMF analysis, or estimate separately
K <- 2  # Adjust based on your population structure analysis

# Run LFMM2 for each environmental variable
for(env_var in 1:ncol(X)) {
  env_name <- colnames(X)[env_var]
  
  # Fit LFMM2 model
  mod_lfmm2 <- lfmm2(input = Y, 
                     env = X[, env_var, drop = FALSE], 
                     K = K,
                     effect.sizes = TRUE)
  
  # Compute p-values and effect sizes
  test_results <- lfmm2.test(object = mod_lfmm2,
                             input = Y,
                             env = X[, env_var, drop = FALSE], 
                             full = TRUE)
  
  pv_lfmm2 <- test_results$pvalues
  effect_sizes <- test_results$effect.sizes
  
  # Store results
  results_list[[env_name]] <- list(
    model = mod_lfmm2,
    pvalues = pv_lfmm2,
    effect_sizes = effect_sizes  # Fixed: use the extracted variable
  )
  
  pvalue_matrices[[env_name]] <- pv_lfmm2
}

# Create Manhattan plots for each environmental variable
par(mfrow = c(ceiling(ncol(X)/2), 2))

for(env_var in names(pvalue_matrices)) {
  pv <- pvalue_matrices[[env_var]]
  
  # Manhattan plot
  plot(-log10(pv), 
       cex = 0.5, 
       xlab = "Locus", 
       ylab = "-Log10(P)", 
       col = "blue",
       main = paste("Manhattan Plot:", env_var))
  
  # Add significance thresholds
  abline(h = -log10(0.05), col = "red", lty = 2, lwd = 1)  # p = 0.05
  abline(h = -log10(0.001), col = "red", lty = 1, lwd = 1)  # p = 0.001
  
  # Add legend
  legend("topright", 
         legend = c("p = 0.05", "p = 0.001"), 
         col = "red", 
         lty = c(2, 1), 
         cex = 0.8)
}

# Reset plotting parameters
par(mfrow = c(1, 1))

# Multiple testing correction and summary

for(env_var in names(pvalue_matrices)) {
  pv <- pvalue_matrices[[env_var]]
  
  # Apply Bonferroni correction
  pv_bonf <- p.adjust(pv, method = "bonferroni")
  
  # Apply FDR correction
  pv_fdr <- p.adjust(pv, method = "fdr")
  
  # Count significant associations
  sig_005 <- sum(pv < 0.05, na.rm = TRUE)
  sig_001 <- sum(pv < 0.001, na.rm = TRUE)
  sig_bonf <- sum(pv_bonf < 0.05, na.rm = TRUE)
  sig_fdr <- sum(pv_fdr < 0.05, na.rm = TRUE)
  
  cat(sprintf("\n%s:\n", env_var))
  cat(sprintf("  Significant at p < 0.05: %d loci (%.2f%%)\n", 
              sig_005, 100 * sig_005 / length(pv)))
  cat(sprintf("  Significant at p < 0.001: %d loci (%.2f%%)\n", 
              sig_001, 100 * sig_001 / length(pv)))
  cat(sprintf("  Significant after Bonferroni: %d loci (%.2f%%)\n", 
              sig_bonf, 100 * sig_bonf / length(pv)))
  cat(sprintf("  Significant after FDR: %d loci (%.2f%%)\n", 
              sig_fdr, 100 * sig_fdr / length(pv)))
  
  # Store corrected p-values
  results_list[[env_var]]$pvalues_bonf <- pv_bonf
  results_list[[env_var]]$pvalues_fdr <- pv_fdr
  
  # Identify top associated loci
  top_loci_idx <- order(pv)[1:min(10, length(pv))]
  top_loci_names <- colnames(Y)[top_loci_idx]
  top_pvalues <- pv[top_loci_idx]
  
  cat("  Top 10 associated loci:\n")
  for(i in 1:length(top_loci_idx)) {
    cat(sprintf("    %s: p = %.2e\n", top_loci_names[i], top_pvalues[i]))
  }
}

# Create combined visualization
# Combine results for multi-environment plot
if(ncol(X) > 1) {
  # Create a comprehensive results dataframe
  all_results <- data.frame()
  
  for(env_var in names(pvalue_matrices)) {
    pv <- pvalue_matrices[[env_var]]
    temp_df <- data.frame(
      Locus = 1:length(pv),
      Locus_name = colnames(Y),
      PValue = pv,
      NegLog10P = -log10(pv),
      Environment = env_var
    )
    all_results <- rbind(all_results, temp_df)
  }
  
  # Multi-environment Manhattan plot
  p_multi <- ggplot(all_results, aes(x = Locus, y = NegLog10P, color = Environment)) +
    geom_point(alpha = 0.6, size = 0.8) +
    facet_wrap(~Environment, scales = "free_y") +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.001), color = "red") +
    labs(x = "Locus Position",
         y = "-Log10(P-value)") +
    theme_minimal() +
    theme(legend.position = "")
  
  print(p_multi)
}

## LFMM2 with entire env matrix ----

# Run LFMM2 with entire environmental matrix

# Fit LFMM2 model with full environmental matrix
mod_lfmm2 <- lfmm2(input = Y, 
                   env = X,  # Use entire environmental matrix
                   K = K,
                   effect.sizes = TRUE)

# Compute p-values and effect sizes for all variables at once
test_results <- lfmm2.test(object = mod_lfmm2,
                           input = Y,
                           env = X, 
                           full = TRUE)

pv_lfmm2 <- test_results$pvalues
effect_sizes <- test_results$effect.sizes

# Manhattan plot showing -log P
plot(-log10(pv_lfmm2), cex = .3, xlab = "Locus",  ylab = "-Log(P)", col = "blue")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")

# FDR control: computing qvalues
qv_lfmm2  <- qvalue::qvalue(pv_lfmm2, fdr.level = 0.2)


# Orange line: Bonferonni multiple testing threshold for significance
abline(h = -log10(0.1/ncol(Y)), col = "orange")

# the most interesting targets
candidates <- which(qv_lfmm2$significant)

# Candidate loci at FDR level = 20%
points(candidates, -log10(pv_lfmm2)[candidates], cex = .9, col = "brown")


# Genomic offset using LEA ----

# Read predicted environmental data
env_pred <- read.csv("predicted_environment.csv",
                     sep = ";",
                     header = TRUE, 
                     row.names = 1)  # Assuming first column has sample names

# Convert to matrix
X_pred <- as.matrix(env_pred)

# Compute GO with all loci
g_offset <- genetic.offset(input = Y, 
                           env = X, 
                           pred.env = X_pred, 
                           scale = TRUE,
                           K = 2)

# Plot genetic offset vs squared  Euclidean  environmental  distance
Delta = X - X_pred
dist_env2 =  rowSums(Delta^2)   
plot(dist_env2, g_offset$offset, 
     xlab ="Squared Euclidean distance",  
     ylab ="Genetic offset", cex = .6, col = "blue")


# Load the relative fitness
rel_fit <- read.csv("relative_fitness.csv",
                     sep = ";",
                     header = TRUE) 

# Plot the relative fitness vs  genetic offset 
plot(g_offset$offset, - rel_fit$x, 
     ylab ="Relative fitness",  
     xlab ="Genetic offset", cex = .6, col = "blue")

# Log the relative fitness
log_rel_fit <- log(rel_fit$x)

# Plot the log relative fitness vs  genetic offset 
plot(g_offset$offset, - log_rel_fit, 
     ylab ="Relative fitness (log)",  
     xlab ="Genetic offset", cex = .6, col = "blue")

cor(g_offset$offset, log_rel_fit)^2


# BayPass ----

# Need g_baypass file in the same directory as the data files
# Make sure g_baypass is executable

./g_baypass -npop 100 -gfile limaxnivalis.baypass.geno -outprefix core_analysis -nthreads 4

./g_baypass -npop 100 -gfile limaxnivalis.baypass.geno -efile limaxnivalis.obsenv.baypass.cov -scalecov -omegafile core_analysis_mat_omega.out -outprefix env_analysis -nthreads 4


#After running the standard model, BayPass produces several output files. The most important for identifying significant SNPs are:
#*_summary_betai_reg.out - Contains Bayes Factors (BF) for each SNP-covariate association
#*_summary_betai.out - Contains detailed association statistics
#*_summary_pi_xtx.out - Contains XtX differentiation statistics

## Identifying Significant SNPs ----

# Using Bayes Factor (dB) Thresholds
# The most common approach is to use a BF threshold of 20 dB
baypass_results <- read.table("env_analysis_summary_betai_reg.out", header = TRUE)

# Set threshold for Bayes Factor
bf_threshold <- 20

# Identify significant SNPs
significant_snps <- baypass_results[baypass_results$BF.dB. > bf_threshold, ]

# View the significant SNPs
head(significant_snps)

# Count how many significant SNPs
nrow(significant_snps)

# Basic Manhattan plot
manhattan_plot <- ggplot(results_with_pos, aes(x = POS, y = DB, color = factor(CHR))) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 15, linetype = "dashed", color = "blue") +
  labs(
    x = "Genomic Position",
    y = "Bayes Factor (dB)",
    title = "Genome-Environment Association Analysis",
    color = "Chromosome"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_grid(. ~ CHR, scales = "free_x", space = "free_x")


# Using XtX values
source("baypass_utils.R")  # BayPass utility functions

# Read XtX values (if using core model results)
xtx <- read.table("env_analysis_summary_pi_xtx.out", header = TRUE)

# Compute 99% threshold for neutral data
threshold_99 <- quantile(xtx$M_XtX, probs = 0.99)
threshold_95 <- quantile(xtx$M_XtX, probs = 0.95)

# Apply this threshold to identify outliers
outlier_snps <- baypass_results[baypass_results$XtX > threshold_99, ]


## Plotting ----

omega=as.matrix(read.table("core_analysis_mat_omega.out"))
plot.omega(omega=omega)

require(corrplot)
cor.mat=cov2cor(omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))

anacore.snp.res=read.table("env_analysis_summary_pi_xtx.out",h=T)
