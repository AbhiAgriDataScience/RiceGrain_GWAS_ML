# ===================  Introduction to Phenotypic Data Processing script for GWAS in Rice ================================================================================================================================================
# This R script is the first step in our GWAS pipeline for rice phenotypic data. 
# The goal of this script is to process and analyze phenotypic data collected from two different locations: Hangzhou (HZ) and Sanya (SY). 
# The key tasks performed in this script include Best Linear Unbiased Prediction (BLUP) estimation, quality control, normality testing, and correlation analysis of traits.
# The BLUP values are computed using a mixed-effects model, where Location is treated as a fixed effect to adjust for environmental influence, and LINE (genotype ID) is treated as a random effect to estimate genetic variance. 
# After BLUP estimation, we compare our computed BLUPs with provided values, assess data distribution, and apply Inverse Normal Transformation (INT) where necessary.

# ============================================== Workflow of this Script ===================================================================
# Step 1: Load Required Packages & Data
# ---------------------------------------------------------------------------------------------------------------------------------------------
#   - Load necessary R packages (lme4, data.table, ggplot2, etc.).
#   - Read the raw phenotypic dataset (GSTP007.pheno).
#   - Identify traits corresponding to Hangzhou (HZ_) and Sanya (SY_).
# ---------------------------------------------------------------------------------------------------------------------------------------------
# Step 2: Compute BLUP Values
# ---------------------------------------------------------------------------------------------------------------------------------------------
#   - Extract common traits present in both locations.
#   - Reshape the data to long format for mixed model fitting.
#   - Fit a mixed model for each trait:
#   TraitValue ∼ Location + (1 ∣ LINE)
#         - Location as a fixed effect
#         - LINE as a random effect
#   -Extract random effect estimates (BLUP values).
#   -Save computed BLUP values to blup_results/computed_BLUP_from_SY_HZ.csv.
# ---------------------------------------------------------------------------------------------------------------------------------------------
# Step 3: Compare Computed BLUP vs Provided BLUP
# ---------------------------------------------------------------------------------------------------------------------------------------------
#   - Load the provided BLUP values.
#   - Merge with computed BLUP values and calculate correlation between them.
#   - Save results to blup_results/blup_correlation_SY_HZ.csv.
# ---------------------------------------------------------------------------------------------------------------------------------------------
# Step 4: Normality Testing
# ---------------------------------------------------------------------------------------------------------------------------------------------
#   - Check if BLUP traits follow a normal distribution using Lilliefors KS test.
#   - Generate histograms and Q-Q plots to visualize distributions.
#   - Save normality results in blup_results/Normality_Test/normality_results_BLUP.csv.
# ---------------------------------------------------------------------------------------------------------------------------------------------
# Step 5: Apply Inverse Normal Transformation (INT)
# ---------------------------------------------------------------------------------------------------------------------------------------------
#   - Apply INT transformation to non-normally distributed traits.
#   - Save transformed dataset as blup_results/pheno_BLUP_INT.csv.
#   - Re-run normality tests and save results in blup_results/Normality_Test/normality_results_INT_BLUP.csv.
# ---------------------------------------------------------------------------------------------------------------------------------------------
# Step 6: Correlation Analysis
# ---------------------------------------------------------------------------------------------------------------------------------------------
#   - Compute the correlation matrix for BLUP traits.
#   - Generate a correlation heatmap.
#   - Identify highly correlated traits and save results to blup_results/High_Correlation_Traits_BLUP.csv.

# ======================================== File Output Summary =================================================================================
# File Name	                                                    Description
# -----------------------------------------------------------   ----------------------------------------------------------------------------------
# blup_results/computed_BLUP_from_SY_HZ.csv	                    Computed BLUP values for each trait.
# blup_results/blup_correlation_SY_HZ.csv	                      Correlation between computed and provided BLUPs.
# blup_results/Normality_Test/normality_results_BLUP.csv	      Normality test results before transformation.
# blup_results/Normality_Test/normality_results_INT_BLUP.csv	  Normality test results after Inverse Normal Transformation (INT).
# blup_results/pheno_BLUP_cleaned.csv	                          Final cleaned BLUP phenotype dataset.
# blup_results/pheno_BLUP_INT.csv	                              Inverse Normal Transformed (INT) BLUP dataset.
# blup_results/Correlation_Matrix_BLUP.csv	                    Trait correlation matrix.
# blup_results/Correlation_Heatmap_BLUP.png	                    Heatmap showing correlation between traits.
# blup_results/High_Correlation_Traits_BLUP.csv	                List of highly correlated traits (> 0.75 correlation).


# -------------------------------------------------------------------------------------------------------------------------------------------------------
# GWAS_Rice - Phenotypic Data Processing
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Author: Abhishek Shrestha
# Date: 5 March 2025
# Description: 
# This script processes rice phenotypic data, computes BLUP values, 
# performs data quality control, checks normality, applies transformations, 
# and analyzes correlations to prepare for GWAS analysis.
# -------------------------------------------------------------------------------------------------------------------------------------------------------

# Step 1: Load required packages
# -------------------------------------------------------------------------------------------------------------------------------------------------------
library(lme4)       # Mixed model fitting
library(data.table) # Faster data handling
library(ggplot2)    # Visualization
library(caret)      # Model evlauation (correlation comparison)
library(reshape2)  # For data reshaping
library(dplyr)
library(tidyr)
library(readr)
library(bestNormalize)
library(nortest)
library(corrplot)

#setwd("~/Documents/R/GWAS_rice") # macbook
setwd("D:/My projects_Oct 2023/GWAS/GWAS_rice") # PC

# Step 2: Compute BLUP Values and Compare with Provided BLUP Data
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Load phenotype data
pheno <- fread("data/GSTP007.pheno", header=TRUE)

# Identify HZ_ (Hangzhou) and SY_ (Sanya) traits
hz_traits <- grep("^HZ_", colnames(pheno), value = TRUE)
sy_traits <- grep("^SY_", colnames(pheno), value = TRUE)

# Ensure the traits match between SY and HZ
common_traits <- gsub("HZ_", "", hz_traits)  # Extract trait names
matching_sy_traits <- paste0("SY_", common_traits)

# Check if all SY_ traits exist in the dataset
existing_sy_traits <- matching_sy_traits %in% colnames(pheno)

# Keep only traits that exist in both SY and HZ
common_traits <- common_traits[existing_sy_traits]
hz_traits <- hz_traits[existing_sy_traits]
matching_sy_traits <- matching_sy_traits[existing_sy_traits]

# Step 3: Compute BLUP Values for Each Trait
# -------------------------------------------------------------------------------------------------------------------------------------------------------

computed_blup <- data.frame(LINE = pheno$LINE) # Create an empty dataframe for computed BLUP values
all_long_data <- data.frame() # Create a master dataframe to store all long-format data

# Loop through each trait and compute BLUP using SY & HZ values
for (i in seq_along(common_traits)) {
  trait_name <- common_traits[i]
  hz_trait <- hz_traits[i]
  sy_trait <- matching_sy_traits[i]
  
  cat("\nComputing BLUP for:", trait_name, "\n")
  
  # Check if both traits exist before proceeding
  if (!(hz_trait %in% colnames(pheno)) | !(sy_trait %in% colnames(pheno))) {
    cat("Skipping", trait_name, "- Missing trait in dataset\n")
    next  # Skip this trait if either HZ_ or SY_ is missing
  }
  
  # Reshape the data to long format for mixed model analysis
  long_data <- melt(pheno, id.vars = "LINE", measure.vars = c(hz_trait, sy_trait),
                    variable.name = "Location", value.name = "TraitValue")
  
  # Fix `Location` column (convert full trait name to just HZ or SY)
  long_data$Location <- ifelse(grepl("^HZ_", long_data$Location), "HZ", "SY")
  
  # Add a column for the trait name
  long_data$Trait <- trait_name  
  
  # Append to master dataframe
  all_long_data <- rbind(all_long_data, long_data)
  
  # Convert `Location` to a factor
  long_data$Location <- as.factor(long_data$Location)
  
  # Remove rows with missing TraitValues
  long_data <- na.omit(long_data)
  
  # Ensure there are at least 3 observations to fit the model
  if (nrow(long_data) < 3) {
    cat("Skipping", trait_name, "- Not enough data for BLUP computation\n")
    next
  }
  
  # Fit a mixed model: TraitValue ~ Location + (1|LINE)
  model <- tryCatch({
    lmer(TraitValue ~ Location + (1|LINE), data = long_data, REML = TRUE) # Location as a fixed effect to remove confounding and LINE as a random effect 
  }, error = function(e) {
    cat("Error in model fitting for", trait_name, "- Skipping this trait\n")
    return(NULL)
  })
  
  # If model fitting failed, skip this trait
  if (is.null(model)) next
  
  
  # Extract BLUP values (random effect estimates for LINE)
  blup_values <- ranef(model)$LINE[,1] 
  observed_mean <- mean(long_data$TraitValue, na.rm=TRUE) # Compute the mean of observed trait values for centering
  
  # Scale BLUP values to match observed trait distribution
  blup_values <- blup_values + observed_mean  # Center around observed mean
  
  # Ensure correct alignment of BLUP values with LINE
  matched_blup <- rep(NA, nrow(pheno))  # Initialize with NAs
  line_names <- rownames(ranef(model)$LINE)  # Get names of included lines
  matched_blup[match(line_names, pheno$LINE)] <- blup_values  # Align to original LINE order
  
  # Store BLUP values in the dataframe
  computed_blup[[paste0("BLUP_", trait_name)]] <- matched_blup
}

# Save computed BLUP values
write.csv(computed_blup, "blup_results/computed_BLUP_from_SY_HZ.csv", row.names=FALSE)
cat("Computed BLUP values saved in 'blup_results/computed_BLUP_from_SY_HZ.csv'\n")


# Step 4: Compare Computed BLUP vs Provided BLUP
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Load provided BLUP data
blup_provided <- fread("data/GSTP007.pheno", select=c("LINE", grep("^BLUP_", colnames(fread("data/GSTP007.pheno")), value = TRUE)))

# Load computed BLUP data
blup_computed <- fread("results/computed_BLUP_from_SY_HZ.csv")

# Merge provided and computed BLUP datasets on LINE
blup_comparison <- merge(blup_provided, blup_computed, by="LINE", all = TRUE)

# Print available columns for debugging
print("Available columns in merged dataset:")
print(colnames(blup_comparison))

# Initialize correlation results dataframe
cor_results <- data.frame(Trait=character(), Correlation=numeric(), stringsAsFactors=FALSE)

# Iterate over each BLUP trait in provided data
for (trait in grep("^BLUP_", colnames(blup_provided), value = TRUE)) {
  # Identify correct column names in merged dataset
  trait_provided <- paste0(trait, ".x")
  trait_computed <- paste0(trait, ".y")
  
  # Check if both columns exist
  if (!(trait_provided %in% colnames(blup_comparison)) || !(trait_computed %in% colnames(blup_comparison))) {
    cat("Skipping", trait, "- No matching computed BLUP found\n")
    next
  }
  
  # Compute correlation
  cor_value <- cor(blup_comparison[[trait_provided]], blup_comparison[[trait_computed]], use="complete.obs")
  
  cor_results <- rbind(cor_results, data.frame(
    Trait = trait,
    Correlation = cor_value
  ))
}

# Save correlation results
write.csv(cor_results, "results/blup_correlation_SY_HZ.csv", row.names=FALSE)


# Step 5: Plot correlation results
# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot correlation results if there is data
if (nrow(cor_results) > 0) {
  p <- ggplot(cor_results, aes(x=Trait, y=Correlation)) +
    geom_bar(stat="identity", fill="blue") +
    coord_flip() +
    ggtitle("Correlation Between Provided and Computed BLUP (SY & HZ)") +
    theme_minimal()
  
  # Save plot
  ggsave("results/blup_correlation_SY_HZ.png", plot = p, width = 10, height = 6)
  
  cat("BLUP correlation analysis completed. Results saved in 'results/blup_correlation_SY_HZ.csv'\n")
} else {
  cat("No valid correlations computed. Check data consistency.\n")
}

# Print results
print(cor_results)

# > print(var(blup_comparison$BLUP_Yield_per_plant.y, na.rm = TRUE))
# [1] 1.598958
# > print(var(blup_comparison$BLUP_Yield_per_plant.x, na.rm = TRUE))
# [1] 1.599012
# Minimal difference in LINE effects

# BLUP Computation is Correct
# 
# Since the variance is nearly identical, your mixed model is producing results very close to the provided values.
# This suggests that:
#   The model structure is fine.
# The BLUPs are capturing the expected genetic variation.
# Minimal Difference in LINE Effects
# 
# Because your random effect variance (LINE) is low, the computed BLUPs stay close to the fixed effect mean.
# However, they still maintain the correct ranking structure.
# Singular Fit (isSingular) in Some Traits
# 
# For some traits, the random effect variance is near 0, causing isSingular warnings.
# This suggests that some traits might not have enough variation to estimate BLUPs accurately.


# Step 6: Assess Fixed and Random Effects Contribution
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# To quantify the effects of Location (fixed effect) and LINE (random effect), you can check:
#   
#   (A) Variance Partitioning (Intraclass Correlation Coefficient - ICC)
# If LINE variance is high, it means genetic differences dominate.
# If Location variance is high, it means the environment dominates.
# If isSingular appears, it suggests LINE variance is near zero.

library(lme4)

trait_results <- data.frame(Trait=character(), Line_Var=numeric(), Residual_Var=numeric(), ICC=numeric(), Location_Effect=numeric(), stringsAsFactors=FALSE)

for (trait in unique(all_long_data$Trait)) {
  subset_data <- subset(all_long_data, Trait == trait)
  
  if (nrow(subset_data) < 3) next  # Skip traits with too few observations
  
  model <- lmer(TraitValue ~ Location + (1 | LINE), data = subset_data, REML = TRUE)
  
  # Extract variance components
  var_comp <- as.data.frame(VarCorr(model))
  line_var <- var_comp$vcov[1]  # LINE variance
  residual_var <- var_comp$vcov[2]  # Residual variance
  total_var <- line_var + residual_var
  
  # Calculate ICC (proportion of variance due to LINE)
  icc <- line_var / total_var
  
  # Extract Location effect size
  loc_effect <- summary(model)$coefficients["LocationSY", "Estimate"]
  
  trait_results <- rbind(trait_results, data.frame(
    Trait = trait,
    Line_Var = line_var,
    Residual_Var = residual_var,
    ICC = icc,
    Location_Effect = loc_effect
  ))
}

# Display results
print(trait_results)

# Save computed BLUP values
write.csv(trait_results, "results/fixed_random_effect.csv", row.names=FALSE)

# 2. Plot influence of LINE vs Location
# Visualize the contribution of LINE vs Residual variance:
library(ggplot2)
ggplot(trait_results, aes(x=Trait, y=ICC)) +
  geom_bar(stat="identity", fill="blue") +
  coord_flip() +
  ggtitle("Proportion of Variation Explained by LINE (ICC)") +
  theme_light()
ggsave("results/ICC_variation_plot.png")

# Plot Location effect across traits:
ggplot(trait_results, aes(x=Trait, y=Location_Effect)) +
  geom_bar(stat="identity", fill="red") +
  coord_flip() +
  ggtitle("Effect of Location on Trait Values") +
  theme_light()
ggsave("results/Location_effect_plot.png")

# 3. Interpret the Results
# If ICC is high (close to 1) → Genetic differences (LINE) explain most of the variation.
# If ICC is low (~0 or isSingular warning) → Little genetic influence, mostly environmental variation.
# If Location_Effect is significant → Location strongly influences trait values.

# Conclusion
# If ICC is high → LINE (genetics) drives variation.
# If Location_Effect is high → Location (environment) drives variation.
# If isSingular occurs → LINE variation is too low, BLUPs are unreliable.


# Trait	              LINE Variance	            Residual Variance	            ICC (Genetic Influence)	              Location Effect
# Yield_per_plant	      High	                    Low	                            Strong LINE effect	                Negligible
# Panicle_number	    Moderate	                Moderate	                        Mixed influence	                    Significant
# Grain_number	      Low	                        High	                          Mostly environment	                    Strong
# Height	              High	                    Low	                                Strong LINE effect	              Negligible



# ----------------------------------------------------------------------------------------------------------------------
# Step 7: Pre-processing & Quality Control of Phenotypic Data
# ----------------------------------------------------------------------------------------------------------------------

# Step 7.1 Load phenotype data (BLUP dataset)
# ----------------------------------------------------------------------------------------------------------------------
pheno_blup <- fread("results/computed_BLUP_from_SY_HZ.csv", header = TRUE)
setDT(pheno_blup)  # Ensure it's a data.table

# Check structure and missing values
str(pheno_blup)
head(pheno_blup)
summary(pheno_blup)
cat("Total missing values:", sum(is.na(pheno_blup)), "\n")
# Total missing values: 8 

# Ensure "LINE" column is present
if (!("LINE" %in% colnames(pheno_blup))) {
  stop("Error: Column 'LINE' not found. Check the dataset column names.")
}

# ----------------------------------------------------------------------------------------------------------------------
#  Step 7.2: Pre-processing & Quality Control of Phenotypic Data
# ----------------------------------------------------------------------------------------------------------------------

# Retain taxa with at least 70% non-missing TraitValue
pheno_filtered <- pheno_blup %>%
  filter(rowMeans(is.na(select(., where(is.numeric)))) <= 0.3)

cat("Remaining unique taxa after filtering:", n_distinct(pheno_filtered$LINE), "\n")
# Remaining unique taxa after filtering: 1494 

# Perform all preprocessing steps
pheno_filtered <- pheno_filtered %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%  # Impute missing values
  mutate(across(where(is.numeric), ~ ifelse(. < 0, min(.[. > 0], na.rm = TRUE) / 2, .))) %>%  # Remove negatives
  mutate(across(where(is.numeric), ~ ifelse(abs(scale(.)) > 3, NA_real_, .))) %>%  # Remove outliers
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))  # Re-impute missing values

# Fix any non-numeric columns
pheno_filtered <- pheno_filtered %>%
  mutate(across(where(is.list), ~ as.numeric(unlist(.)))) %>%
  mutate(across(where(is.matrix), ~ as.numeric(.)))

# Save cleaned dataset
write_csv(pheno_filtered, "blup_results/pheno_BLUP_cleaned.csv")
cat("Final BLUP phenotype dataset saved as 'blup_results/pheno_BLUP_cleaned.csv'\n")

## ----------------------------------------------------------------------------------------------------------------------
# Step 7.3: Normality Test (Check if Traits Are Normally Distributed)
# ----------------------------------------------------------------------------------------------------------------------

# Define output directory for plots
plot_dir <- "blup_results/Normality_Test/Plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Initialize normality results dataframe
normality_results <- data.frame(Trait = character(), P_Value = numeric(), Normal = character(), stringsAsFactors = FALSE)

# Perform normality test for each trait (excluding "LINE")
for (trait in colnames(pheno_filtered)[-1]) {  
  data <- pheno_filtered[[trait]]
  data <- data[is.finite(data)]  # Remove non-finite values
  
  # Perform Lilliefors KS normality test
  ks_test <- lillie.test(data)
  
  # Store results
  normality_results <- rbind(normality_results, 
                             data.frame(Trait = trait, P_Value = ks_test$p.value, 
                                        Normal = ifelse(ks_test$p.value > 0.05, "Yes", "No")))
  
  # Generate histogram
  hist_plot <- ggplot(pheno_filtered, aes(x = .data[[trait]])) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = paste("Histogram of", trait), x = "Trait Value", y = "Frequency") +
    theme_light()
  
  ggsave(filename = paste0(plot_dir, "/histogram_", trait, ".png"), plot = hist_plot, width = 6, height = 4, dpi = 300)
  
  # Generate Q-Q Plot
  png(paste0(plot_dir, "/qqplot_", trait, ".png"), width = 800, height = 600)
  qqnorm(data, main = paste("Q-Q Plot of", trait))
  qqline(data, col = "red")
  dev.off()
}

# Save normality test results
write_csv(normality_results, "blup_results/Normality_Test/normality_results_BLUP.csv")

cat("Normality check completed. Review histograms and Q-Q plots to decide transformations.\n")

# only BLUP_Panicle_length was found to be normally distributed (0.055043874285982924,Yes)

# ----------------------------------------------------------------------------------------------------------------------
# Step 7. 4: Apply Inverse Normal Transformation (INT) for Non-Normal Traits
# ----------------------------------------------------------------------------------------------------------------------

# Select traits that failed normality test
non_normal_traits <- normality_results$Trait[normality_results$Normal == "No"]

# Apply INT transformation to non-normal traits
pheno_int <- pheno_filtered

for (trait in non_normal_traits) {
  pheno_int[[trait]] <- qnorm((rank(pheno_filtered[[trait]], ties.method = "average") - 0.5) / length(pheno_filtered[[trait]]))
}

# Save transformed dataset
write_csv(pheno_int, "blup_results/pheno_BLUP_INT.csv")
cat("Inverse Normal Transformation applied. Saved as 'blup_results/pheno_BLUP_INT.csv'\n")

# Perform Normality test (After INT transformation)

# Initialize normality results dataframe (Post-INT)
normality_results_int <- data.frame(Trait = character(), P_Value = numeric(), Normal = character(), stringsAsFactors = FALSE)

# Perform normality test for transformed traits
for (trait in non_normal_traits) {  
  data <- pheno_int[[trait]]
  data <- data[is.finite(data)]  # Remove non-finite values
  
  # Perform Lilliefors KS normality test
  ks_test <- lillie.test(data)
  
  # Store results
  normality_results_int <- rbind(normality_results_int, 
                                 data.frame(Trait = trait, P_Value = ks_test$p.value, 
                                            Normal = ifelse(ks_test$p.value > 0.05, "Yes", "No")))
  
  # Generate histogram
  hist_plot <- ggplot(pheno_int, aes(x = .data[[trait]])) +
    geom_histogram(bins = 30, fill = "lightcoral", color = "black", alpha = 0.7) +
    labs(title = paste("Histogram of", trait, "(After INT)"), x = "Trait Value", y = "Frequency") +
    theme_light()
  
  ggsave(filename = paste0(plot_dir, "/histogram_INT_", trait, ".png"), plot = hist_plot, width = 6, height = 4, dpi = 300)
  
  # Generate Q-Q Plot
  png(paste0(plot_dir, "/qqplot_INT_", trait, ".png"), width = 800, height = 600)
  qqnorm(data, main = paste("Q-Q Plot of", trait, "(After INT)"))
  qqline(data, col = "red")
  dev.off()
}

# Save new normality results
write_csv(normality_results_int, "blup_results/Normality_Test/normality_results_INT_BLUP.csv")

cat("Post-INT Normality check completed. Review histograms and Q-Q plots to confirm transformation effectiveness.\n")
# Except for BLUP_Panicle_number (5.510536666596535e-10,No), all other traits are normally distributed now after INT transformation

# ----------------------------------------------------------------------------------------------------------------------
# Step 7.5: Correlation Analysis (BLUP Traits)
# ----------------------------------------------------------------------------------------------------------------------

# Compute correlation matrix excluding LINE column
cor_matrix <- cor(pheno_int %>% select(-LINE), use = "pairwise.complete.obs")

# Save correlation matrix
write_csv(as.data.frame(cor_matrix), "blup_results/Correlation_Matrix_BLUP.csv")

# Generate and save correlation heatmap
png("blup_results/Correlation_Heatmap_BLUP.png", width = 800, height = 600)
corrplot(cor_matrix, method = "color", type = "lower", tl.cex = 0.8, col = colorRampPalette(c("blue", "white", "red"))(200))
dev.off()

# Identify highly correlated trait pairs
high_corr_traits <- which(abs(cor_matrix) > 0.75 & abs(cor_matrix) < 1, arr.ind = TRUE)

# Extract trait pairs and their correlation values
high_corr_traits_df <- data.frame(
  Trait1 = rownames(cor_matrix)[high_corr_traits[, 1]],
  Trait2 = colnames(cor_matrix)[high_corr_traits[, 2]],
  Correlation = cor_matrix[high_corr_traits]
) %>%
  distinct()

# Save highly correlated traits
write_csv(high_corr_traits_df, "blup_results/High_Correlation_Traits_BLUP.csv")
cat("Correlation analysis completed. High correlation traits saved.\n")

# We do not have any traits that are highly correlated (>0.75).

cat("Phenotypic data processing completed. Outputs are saved in 'blup_results'.\n")

# ----------------------------------------------------------------------------------------------------------------------
# Pre-processing processed phenotype data

pheno_INT_data = fread("blup_results/pheno_BLUP_INT.csv", header = TRUE)
str(pheno_INT_data)
# except for panicle length, all traits (9) are INT transformed 

# Rename columns by removing "BLUP_" prefix
colnames(pheno_INT_data) = gsub("^BLUP_", "", colnames(pheno_INT_data))

# Remove individuals with missing phenotype data
pheno_INT_data <- na.omit(pheno_INT_data)

cat("Phenotype data dimensions:", dim(pheno_INT_data), "\n")

nrow(pheno_INT_data)

# Fix column name in pheno_data to match GAPIT's format
colnames(pheno_INT_data)[1] = "Taxa"

str(pheno_INT_data)
head(pheno_INT_data)

sum(is.na(pheno_INT_data))

# Save the processed phenotype data (trait name corrected and free of NA values) for GAPIT analysis
fwrite(pheno_INT_data, "GAPIT_data/pheno_BLUP_INT_cleaned.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Processed phenotype data saved as 'GAPIT_data/pheno_BLUP_INT_cleaned.txt'\n")

write.csv(pheno_INT_data, "GAPIT_data/pheno_BLUP_INT_cleaned.csv", row.names = FALSE)
