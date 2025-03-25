#!/bin/bash

# Define input and output filenames
GWAS_RESULTS="pGWAS/Grain_weight_GWAS_prepared_for_PLINK.txt"
PLINK_BFILE="data/GSTP007.bed/1495Hybrid_MSUv7"  # Change this to your PLINK binary file prefix
OUT_FILE="pGWAS/lead_SNPs/lead_SNPs"

# Run PLINK clumping
plink --bfile $PLINK_BFILE \
      --clump $GWAS_RESULTS \
      --clump-p1 0.001 \
      --clump-r2 0.1 \ 
      --clump-kb 250 \ 
      --out $OUT_FILE

# Check if clumping was successful
if [ -f "${OUT_FILE}.clumped" ]; then
    echo "SNP Clumping Completed! Output saved in ${OUT_FILE}.clumped"
else
    echo "SNP Clumping Failed. Check your PLINK input files!"
fi
