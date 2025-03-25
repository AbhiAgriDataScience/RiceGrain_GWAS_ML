library(gplots)         # For heatmaps and plots
library(scatterplot3d)  # For 3D visualization
library(multtest)
library(MASS)
library(GAPIT)
library(rrBLUP) # Required for kinship calculation

# Load libraries
library(data.table)
library(ggplot2)
library(qqman)
library(GAPIT)

library(tidyr)
library(dplyr)
library(readr)



library(nortest)   # Lilliefors KS test for normality

library(ggcorrplot)
library(corrplot)


library(snpStats)       # SNP genotype/PLINK file processing
library(data.table)     # Efficient data loading
library(dplyr)          # Data manipulation
library(tidyverse)      # Data manipulation
library(HardyWeinberg)  # HWE Test
library(tidyr)          # Reshaping data
library(readxl)         # Reading Excel files

# ===============

library(GAPIT)

# Set working directory
setwd("~/Documents/R/GWAS_rice") # Macbook
setwd("D:/My projects_Oct 2023/GWAS/GWAS_rice")

# ----------------------------------------------------------------------------------------------------------------------
# # Step 1: Load Phenotype Data
# pheno_data <- read.csv("genotype_qc_analysis/filtered_phenotype.csv", header=TRUE)
# Step 1: Load and prepare phenotype data
pheno_data = read.csv("GAPIT_data/pheno_BLUP_INT_cleaned.csv", header = TRUE)

head(pheno_data)
cat("Phenotype data dimensions:", dim(pheno_data), "\n")
nrow(pheno_data)
str(pheno_data)


# ----------------------------------------------------------------------------------------------------------------------
# To cross check just to make sure that individual IDs in both phenotype data and genotype data are identical
# IDENTIFY INDIVIDUAKL IDS IN BOTH DATASETS

# Read the genotype HapMap file (without modifying headers)
geno_hmp <- fread("genotype_qc_analysis/gwas_ready.hmp.txt", header = TRUE, data.table = FALSE)

print(head(geno_hmp)[, 1: 50]) # Display the first few rows and few columnsto check the structure

# Extract inbdividual IDs 
# Extract individual IDs from genotype data (column names after metadata columns)
geno_individuals <- colnames(geno_hmp)[12:ncol(geno_hmp)]

# Extract individual IDs from phenotype data (assuming the first column is the individual ID)
pheno_individuals <- pheno_data[, 1]  # Modify this if the column has a different name

# Check their lengths
cat("Number of individuals in genotype data:", length(geno_individuals), "\n")
cat("Number of individuals in phenotype data:", length(pheno_individuals), "\n")

# Identifuy missing individuals

# Identify individuals in genotype but not in phenotype
missing_in_pheno <- setdiff(geno_individuals, pheno_individuals)

# Identify individuals in phenotype but not in genotype
missing_in_geno <- setdiff(pheno_individuals, geno_individuals)

# Print results
cat("Individuals present in genotype but missing in phenotype:", missing_in_pheno, "\n")
cat("Individuals present in phenotype but missing in genotype:", missing_in_geno, "\n")

# The genotype file has the same indivudal IDs as phentoype data. Therefore, we can proceed with the GAPIT analysis.

# ----------------------------------------------------------------------------------------------------------------------

# Read the Hapmap file
geno_hmp <- fread("genotype_qc_analysis/gwas_ready.hmp.txt", header = FALSE, data.table = FALSE)

# Display the first few rows
print(head(geno_hmp)[, 1:50])

dim(geno_hmp)
#[1] 484881   1505
# ----------------------------------------------------------------------------------------------------------------------
# GAPIT analysis (HOPEFULLY THIS WORKS)
library(GAPIT)

myGAPIT <- GAPIT(
  Y = pheno_data[, 1:2],     # Phenotype data
  G = geno_hmp,       # Genotype data
  PCA.total = 3,      # Use first 3 principal components
  model = "MLM",      # Mixed Linear Model (MLM) for population structure
  SNP.MAF = 0.05,
  file.output = TRUE  # Save results
)

# GWAS for Panicle_length
Panicle_length_GWAS <- GAPIT(
  Y = pheno_data[, c("Taxa", "Panicle_length")],     # Phenotype data
  G = geno_hmp,       # Genotype data
  PCA.total = 3,      # Use first 3 principal components
  model = "MLM",      # Mixed Linear Model (MLM) for population structure
  SNP.MAF = 0.05,
  file.output = TRUE  # Save results
)

# GWAS for Grain_number_per_panicle
Grain_number_per_panicle_GWAS <- GAPIT(
  Y = pheno_data[, c("Taxa", "Grain_number_per_panicle")],     # Phenotype data
  G = geno_hmp,       # Genotype data
  PCA.total = 3,      # Use first 3 principal components
  model = "MLM",      # Mixed Linear Model (MLM) for population structure
  SNP.MAF = 0.05,
  file.output = TRUE  # Save results
)

# GWAS for Height
Height_GWAS <- GAPIT(
  Y = pheno_data[, c("Taxa", "Height")],     # Phenotype data
  G = geno_hmp,       # Genotype data
  PCA.total = 3,      # Use first 3 principal components
  model = "MLM",      # Mixed Linear Model (MLM) for population structure
  SNP.MAF = 0.05,
  file.output = TRUE  # Save results
)

# GWAS for Panicle_number
Panicle_number_GWAS <- GAPIT(
  Y = pheno_data[, c("Taxa", "Panicle_number")],     # Phenotype data
  G = geno_hmp,       # Genotype data
  PCA.total = 3,      # Use first 3 principal components
  model = "MLM",      # Mixed Linear Model (MLM) for population structure
  SNP.MAF = 0.05,
  file.output = TRUE  # Save results
)
# There were none SNPs for this trait
