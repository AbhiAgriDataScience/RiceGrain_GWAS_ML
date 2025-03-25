#!/bin/bash

set -e  # Exit immediately if a command fails

# Define input and output file paths
VCF_INPUT="genotype_qc_analysis/final_genotype_fixed.vcf"  # Your fixed VCF file
HMP_OUTPUT="genotype_qc_analysis/gwas_ready.hmp.txt"  # Output Hapmap file

echo "Starting VCF to Hapmap conversion..."

# Check if TASSEL is installed in the Conda environment
if ! command -v run_pipeline.pl &> /dev/null; then
    echo "Error: TASSEL is not found in your Conda environment."
    echo "Try activating the Conda environment with: conda activate plink_env"
    exit 1
fi

# Run TASSEL pipeline for conversion
run_pipeline.pl -Xmx8g \
    -vcf ${VCF_INPUT} \
    -export ${HMP_OUTPUT%.hmp.txt} \
    -exportType Hapmap

echo "Hapmap file created: ${HMP_OUTPUT}"
