#!/bin/bash

# ------------------------------
# Step 1: Define Input & Output
# ------------------------------
LEAD_SNPS_FILE="pGWAS/lead_SNPs/lead_SNPs.clumped"  # Output from PLINK clumping
PLINK_BFILE="data/GSTP007.bed/1495Hybrid_MSUv7"     # Your PLINK binary file
OUT_DIR="pGWAS/finemap_results"                     # Output directory for FINEMAP
N_SAMPLES=$(wc -l data/GSTP007.bed/1495Hybrid_MSUv7.fam | awk '{print $1}') # Count samples in PLINK dataset

# Create output directory if it doesn't exist
mkdir -p $OUT_DIR

# ------------------------------
# Step 2: Extract Lead SNPs for Fine-Mapping
# ------------------------------
echo "Extracting lead SNPs..."
awk 'NR>1 {print $3}' $LEAD_SNPS_FILE > $OUT_DIR/lead_snps.txt
if [ $? -ne 0 ]; then
    echo "Error: Could not extract lead SNPs from $LEAD_SNPS_FILE"
    exit 1
fi
echo "Lead SNP extraction complete!"

# ------------------------------
# Step 3: Generate LD Matrix for Fine-Mapping
# ------------------------------
echo "Generating LD matrix..."
plink --bfile $PLINK_BFILE \
      --r2 square gz \
      --extract $OUT_DIR/lead_snps.txt \
      --out $OUT_DIR/LD_matrix

if [ $? -ne 0 ]; then
    echo "Error: LD matrix generation failed!"
    exit 1
fi
echo "LD matrix generated!"

# ------------------------------
# Step 4: Add Headers to LD Matrix File
# ------------------------------
echo "Adding headers to LD matrix..."
gunzip -c $OUT_DIR/LD_matrix.ld.gz > $OUT_DIR/LD_matrix.ld  # Unzip LD matrix file

# Read SNP IDs
SNP_IDS=($(awk '{print $1}' $OUT_DIR/lead_snps.txt))

# Create the header row (column names)
HEADER="SNP"
for SNP in "${SNP_IDS[@]}"; do
    HEADER="$HEADER\t$SNP"
done
echo -e "$HEADER" > $OUT_DIR/LD_matrix_with_header.ld

# Add SNP row names & append LD matrix data
awk -v snps="${SNP_IDS[*]}" '
BEGIN { split(snps, snp_array, " ") }
{
    print snp_array[FNR] "\t" $0
}' $OUT_DIR/LD_matrix.ld >> $OUT_DIR/LD_matrix_with_header.ld

# Compress back to .gz format
gzip -f $OUT_DIR/LD_matrix_with_header.ld
mv $OUT_DIR/LD_matrix_with_header.ld.gz $OUT_DIR/LD_matrix.ld.gz

if [ $? -ne 0 ]; then
    echo "Error: Failed to add header to LD matrix!"
    exit 1
fi
echo "LD matrix header added successfully!"

# ------------------------------
# Step 5: Prepare FINEMAP Input Files
# ------------------------------
echo "Preparing FINEMAP input files..."

# Create SNP file (SNP information) with headers
echo -e "SNP\tCHR\tBP\tP\tN" > $OUT_DIR/finemap_snps.txt
awk 'NR>1 {print $3, $1, $4, $5, $6}' $LEAD_SNPS_FILE >> $OUT_DIR/finemap_snps.txt

if [ $? -ne 0 ]; then
    echo "Error: Failed to create finemap_snps.txt"
    exit 1
fi

# Create Z-score file (GWAS summary statistics) with headers
echo -e "SNP\tP\tZ" > $OUT_DIR/finemap_z.txt
awk 'NR>1 {print $3, $5, $7}' $LEAD_SNPS_FILE >> $OUT_DIR/finemap_z.txt

if [ $? -ne 0 ]; then
    echo "Error: Failed to create finemap_z.txt"
    exit 1
fi

# Create configuration file for FINEMAP
cat <<EOF > $OUT_DIR/finemap_config.txt
z;${OUT_DIR}/finemap_z.txt
ld;${OUT_DIR}/LD_matrix.ld.gz
snp;${OUT_DIR}/finemap_snps.txt
n;${N_SAMPLES}
EOF

if [ $? -ne 0 ]; then
    echo "Error: Failed to create FINEMAP configuration file!"
    exit 1
fi

echo "FINEMAP input files prepared successfully!"
