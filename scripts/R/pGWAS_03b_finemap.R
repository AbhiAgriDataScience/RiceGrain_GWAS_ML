# This script is designed for Bayesian fine-mapping of genetic associations using SuSiE (Sum of Single Effects model). 
# It processes GWAS results and linkage disequilibrium (LD) matrices to identify credible SNP sets associated with a trait.
# The goal is to prioritize SNPs that are most likely to have a causal effect on the trait, reducing false positives in genome-wide association studies (GWAS).

# Load required libraries
library(data.table)
library(dplyr)
library(susieR)
library(R.utils)

# ------------------------------
# Step 1: Load Input Files
# ------------------------------
lead_snps_file <- "pGWAS/lead_SNPs/lead_SNPs.clumped"
ld_matrix_file <- "pGWAS/finemap_results/LD_matrix.ld.gz"

# Read lead SNPs from PLINK output
lead_snps <- fread(lead_snps_file, header = TRUE)

# Read LD matrix (assumed to be in PLINK format)
ld_matrix <- fread(ld_matrix_file)

# ------------------------------
# Step 2: Prepare FINEMAP Input Data
# ------------------------------
# Use data.table syntax 
finemap_data <- lead_snps[, .(SNP, CHR, BP, P, TOTAL)]  # Select required columns

# Rename columns in data.table
setnames(finemap_data, c("SNP", "CHR", "BP", "P", "TOTAL"), c("SNP", "CHR", "BP", "P", "N"))

# Estimate standard error
finemap_data[, SE := sqrt(1/N)]

# Check first few rows
print(head(finemap_data))

# ------------------------------
# Step 3: Run Fine-Mapping using SuSiE
# ------------------------------
# Extract SNP names (first column)
snp_names <- ld_matrix[[1]]  # First column contains SNP names

# Remove the first column to get only numeric values
X_numeric <- as.matrix(ld_matrix[, -1, with = FALSE])  

# Convert all elements to numeric
X_numeric <- apply(X_numeric, 2, as.numeric)  

# Ensure dimensions match before assigning row/column names
if (length(snp_names) == nrow(X_numeric) && length(snp_names) == ncol(X_numeric)) {
  rownames(X_numeric) <- snp_names  # Assign SNP names as row names
  colnames(X_numeric) <- snp_names  # Assign SNP names as column names
} else {
  stop("Mismatch between SNP names and LD matrix dimensions!")
}

# Convert GWAS p-values to numeric (ensure it's a vector)
y_numeric <- as.numeric(-log10(finemap_data$P))

# Run SuSiE fine-mapping
susie_result <- susie(X_numeric, y_numeric, L = 20) # L = 20 instead of 10 (default)

# Print the summary of the results
summary(susie_result)


# ------------------------------
# Step 4: Extract Credible SNP Sets
# ------------------------------
credible_snps <- susie_result$sets$cs
credible_snps_df <- data.frame(SNP = credible_snps, PIP = susie_result$pip)

# ------------------------------
# Step 5: Save Output
# ------------------------------
output_file <- "pGWAS/finemap_results/finemap_output.csv"
fwrite(credible_snps_df, output_file, sep = "\t", row.names = FALSE)

cat("Bayesian fine-mapping complete! Results saved in:", output_file)


