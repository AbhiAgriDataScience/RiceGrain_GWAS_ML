#!/bin/bash

set -e  # Exit immediately if a command fails

# Define working directory
ANALYSIS_DIR="genotype_qc_analysis"
mkdir -p ${ANALYSIS_DIR}  # Create directory if it doesn't exist

# Input and output file paths
GENOTYPE_PREFIX="data/GSTP007.bed/1495Hybrid_MSUv7"  # Change this to your actual PLINK binary prefix
PHENOTYPE_FILE="blup_results/computed_BLUP_from_SY_HZ.csv"
QC_GENOTYPE_PREFIX="${ANALYSIS_DIR}/qc_filtered_1495Hybrid_MSUv7"
FINAL_GENOTYPE_PREFIX="${ANALYSIS_DIR}/final_genotype"
VCF_OUTPUT="${ANALYSIS_DIR}/final_genotype.vcf"
FIXED_VCF_OUTPUT="${ANALYSIS_DIR}/final_genotype_fixed.vcf"
INDIVIDUALS_FILE="${ANALYSIS_DIR}/filtered_individuals.txt"
SAMPLE_MAP_FILE="${ANALYSIS_DIR}/sample_name_map.txt"

echo "Starting Genotype QC Analysis in ${ANALYSIS_DIR}..."

# Step 1: Perform Genotype Quality Control (QC)
echo "Performing Genotype QC..."
plink --bfile ${GENOTYPE_PREFIX} \
      --geno 0.1 \
      --maf 0.05 \
      --hwe 1e-6 \
      --make-bed \
      --out ${QC_GENOTYPE_PREFIX}

echo "Genotype QC completed. Output saved in ${ANALYSIS_DIR}"

# Step 2: Extract Sample IDs from Phenotype Data (Ensure Proper PLINK Format)
echo "Extracting individuals with phenotype data..."
awk -F',' 'NR>1 {gsub(/"/, "", $1); print $1, $1}' ${PHENOTYPE_FILE} > ${INDIVIDUALS_FILE}

# Remove any trailing spaces, extra characters, and ensure Unix line endings
sed -i.bak 's/"//g' ${INDIVIDUALS_FILE}  # Remove double quotes
dos2unix ${INDIVIDUALS_FILE} 2>/dev/null || true  # Ensure Unix format (ignore errors)
awk '{gsub(/^[ \t]+|[ \t]+$/, ""); print $1, $2}' ${INDIVIDUALS_FILE} > temp.txt
mv temp.txt ${INDIVIDUALS_FILE}

echo "Filtered individuals list saved to ${INDIVIDUALS_FILE}"

# Step 3: Filter Genotype Data to Keep Only Individuals with Phenotypes
echo "Filtering genotype data to match phenotype individuals..."
plink --bfile ${QC_GENOTYPE_PREFIX} \
      --keep ${INDIVIDUALS_FILE} \
      --make-bed \
      --out ${FINAL_GENOTYPE_PREFIX}

echo "Filtered genotype data saved in ${ANALYSIS_DIR}"

# Step 4: Convert Cleaned Genotype Data to VCF Format
echo "Converting cleaned genotype data to VCF..."
plink --bfile ${FINAL_GENOTYPE_PREFIX} \
      --recode vcf \
      --out ${FINAL_GENOTYPE_PREFIX}

echo "VCF file saved as ${VCF_OUTPUT}"

# Step 5: Fix Sample Names in VCF (Remove '_FID' format like Z1_Z1 → Z1)
echo "Fixing sample names in VCF..."
bcftools query -l ${VCF_OUTPUT} | sed 's/\(.*\)_\(.*\)/\1/' > ${SAMPLE_MAP_FILE}

# Verify the mapping file
echo "Sample name mapping file created:"
head ${SAMPLE_MAP_FILE}

# Apply the mapping to correct the VCF header
bcftools reheader -s ${SAMPLE_MAP_FILE} -o ${FIXED_VCF_OUTPUT} ${VCF_OUTPUT}

echo "Fixed VCF file saved as ${FIXED_VCF_OUTPUT}"

echo "Pipeline completed successfully! All results are in ${ANALYSIS_DIR}"
