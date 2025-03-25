# Step 1: Extract Effect and Reference Alleles & Merge with GWAS Results



setwd("D:/My projects_Oct 2023/GWAS/GWAS_rice")

# Load required libraries
library(data.table)
library(dplyr)
library(stringr)

# Define file paths
genotype_file <- "genotype_qc_analysis/gwas_ready.hmp.txt"  # genotype hapmap data file
gwas_file <- "GWAS_GAPIT/GWAS_individual_results/Grain_weight_GWAS_results.csv"    # GWAS result from GAPIT of Grain weight trait

# Load Genotype Data (HapMap format)
# Read the genotype HapMap file (without modifying headers)
genotype_data <- fread(genotype_file, header = TRUE, data.table = FALSE)
print(head(genotype_data)[, 1: 50]) # Display the first few rows and few columnsto check the structure

# Load GWAS Results
gwas_results <- fread(gwas_file)
head(gwas_results)

# Ensure column names are consistent
colnames(genotype_data)[1] <- "SNP"  # Standardize SNP column name
colnames(gwas_results)[2] <- "SNP"   # Ensure SNP column name is the same in both files

# Extract effect/reference alleles from genotype HapMap data
extract_alleles <- function(geno_df) {
  geno_df <- geno_df %>% mutate(
    Reference_Allele = sapply(strsplit(as.character(geno_df$alleles), "/"), `[`, 1),  # First allele
    Effect_Allele = sapply(strsplit(as.character(geno_df$alleles), "/"), `[`, 2)      # Second allele
  )
  
  return(geno_df[, c("SNP", "Reference_Allele", "Effect_Allele")])  # Keep only relevant columns
}

# Apply extraction function
alleles_info <- extract_alleles(genotype_data)

# Merge GWAS results with allele information
gwas_with_alleles <- merge(gwas_results, alleles_info, by = "SNP", all.x = TRUE)
colnames(gwas_with_alleles)

# Estimate missing values
gwas_with_alleles <- gwas_with_alleles %>%
  mutate(
    Beta = -log10(P.value) / (1 - 2 * MAF),  # Estimate Beta
    SE = Beta / -log10(P.value)  # Estimate Standard Error
  )

head(gwas_with_alleles)

# Remove V1 column (rowname turned into column)
# gwas_with_alleles = gwas_with_alleles[, -c("V1")]

# Select and rename columns using data.table syntax
gwas_with_alleles <- gwas_with_alleles[, .(
  CHR = Chr,            # Rename Chr to CHR
  BP = Pos,             # Rename Pos to BP
  SNP = SNP,            # Keep SNP column
  P = P.value,          # Rename P.value to P
  A1 = Effect_Allele,   # Rename Effect_Allele to A1 (effect allele)
  A2 = Reference_Allele,# Rename Reference_Allele to A2 (reference allele)
  Beta = Beta,          # Keep Beta (effect size)
  SE = SE,              # Keep SE (standard error)
  MAF = MAF             # Keep MAF (minor allele frequency)
)]

# Save merged file
fwrite(gwas_with_alleles, "pGWAS/Grain_weight_GWAS_prepared_for_PLINK.csv", sep = ",", row.names = FALSE)
fwrite(gwas_with_alleles, "pGWAS/Grain_weight_GWAS_prepared_for_PLINK.txt", sep = "\t", row.names = FALSE, quote = FALSE)
print("✅ GWAS results merged with effect and reference alleles. File saved as GWAS_results_with_alleles.csv")
