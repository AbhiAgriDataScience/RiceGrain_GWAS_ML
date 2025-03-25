#!/bin/bash

set -e  # Exit on error

# Define file paths
VCF_INPUT="genotype_qc_analysis/final_genotype.vcf"
VCF_OUTPUT="genotype_qc_analysis/final_genotype_fixed.vcf"
SAMPLE_MAP="genotype_qc_analysis/sample_name_map.txt"

echo "Fixing sample names in VCF..."

# Step 1: Generate the mapping file (old name → new name)
bcftools query -l ${VCF_INPUT} | sed 's/\(.*\)_\(.*\)/\1/' > ${SAMPLE_MAP}

# Step 2: Verify mapping file
echo "Sample name mapping file created:"
head ${SAMPLE_MAP}

# Step 3: Apply the new sample names using bcftools
bcftools reheader -s ${SAMPLE_MAP} -o ${VCF_OUTPUT} ${VCF_INPUT}

echo "Fixed VCF file saved as ${VCF_OUTPUT}"
