# -----------------------------------------------------------------------------------------------------
# Step 1: Set Working Directory & Load Required Libraries
# -----------------------------------------------------------------------------------------------------
setwd("D:/My projects_Oct 2023/GWAS/GWAS_rice")
BiocManager::install("STRINGdb", ask= FALSE, force = TRUE)
# Load libraries
library(biomaRt)
library(clusterProfiler)
library(org.Osativa.eg.db)  # Change if using a different species
library(data.table)
library(dplyr)
library(ggplot2)
library(KEGGREST)
library(GenomicRanges)
library(snpStats)  # For SNP annotation
library(STRINGdb) # FOr PPI network analysis

# -----------------------------------------------------------------------------------------------------
# Step 2: Define Input & Output Files
# -----------------------------------------------------------------------------------------------------
finemap_snps_file <- "pGWAS/finemap_results/finemap_snps.txt"
finemap_pip_file <- "pGWAS/finemap_results/finemap_output.csv"
output_dir <- "pGWAS/candidate_genes_new"
genome_db <- "osativa_eg_gene"  # Change for species
window <- 50000  # Search ±50kb for nearby genes
extended_window <- 100000  # Search ±100kb for nearby genes

# Create output directory if it does not exist
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursiv=TRUE)
}

# -----------------------------------------------------------------------------------------------------
# Step 3: Load and Process SNP Data
# -----------------------------------------------------------------------------------------------------
# Step 3.1: Load SNP data
# --------------------------------------------------------------------------------------------------------
finemap_snps <- fread(finemap_snps_file, sep = ' ', header = FALSE)  # SNP, CHR, BP, P, N
setnames(finemap_snps, c("SNP", "CHR", "BP", "P", "N"))  # Rename columns

# --------------------------------------------------------------------------------------------------------
# Step 3.2: Load PIP data
finemap_pip <- fread(finemap_pip_file)  # Columns: SNP, PIP
if (nrow(finemap_snps) != nrow(finemap_pip)) {
  stop("ERROR: The number of SNPs and PIP values do not match!")
}

# Merge SNPs with PIP
finemap_snps[, PIP := finemap_pip$PIP]
fwrite(finemap_snps, file.path(output_dir, "finemap_results_with_PIP.csv"), sep = "\t")

# -----------------------------------------------------------------------------------------------------
# Step 4: Extract Credible SNPs (PIP > 0.95)
# -----------------------------------------------------------------------------------------------------
# Filter SNPs with high posterior probability (PP > 0.95)
credible_snps <- finemap_snps[PIP > 0.95, .(SNP, CHR, BP, PIP)]

# Save credible SNPs
fwrite(credible_snps, file.path(output_dir, "credible_snps.csv"), sep = "\t")

# -----------------------------------------------------------------------------------------------------
# Step 5: Retrieve Gene Annotations
# -----------------------------------------------------------------------------------------------------
# Connect to Ensembl Plants
mart <- useMart("plants_mart", dataset = genome_db, host = "https://plants.ensembl.org")
listDatasets(mart)

# --------------------------------------------------------------------------------------------------------
# Step 5.1: Fetch genes within ±50kb of finemap SNPs
genes <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "external_gene_name"),
  filters = c("chromosome_name", "start", "end"),
  values = list(as.character(finemap_snps$CHR),
                as.numeric(finemap_snps$BP - window),
                as.numeric(finemap_snps$BP + window)),
  mart = mart
)

# Save candidate genes
fwrite(genes, file.path(output_dir, "finemap_candidate_genes.csv"), sep = "\t")
print("Candidate genes saved!")

# --------------------------------------------------------------------------------------------------------
# Step 5.2: Fetch genes within ±50kb of finemap SNPs
genes <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "external_gene_name"),
  filters = c("chromosome_name", "start", "end"),
  values = list(as.character(credible_snps$CHR),
                as.numeric(credible_snps$BP - window),
                as.numeric(credible_snps$BP + window)),
  mart = mart
)

# Save candidate genes
fwrite(genes, file.path(output_dir, "credibleSNPs_candidate_genes.csv"), sep = "\t")
print("Candidate genes saved!")

# --------------------------------------------------------------------------------------------------------
# Step 6: SNP Functional Annotation
# --------------------------------------------------------------------------------------------------------
# Step 6.A: Finemap SNPS
# Step 6.A.1: Define GRanges object for SNPs (all fine-mapped SNPs)
snp_ranges <- GRanges(
  seqnames = paste0("Chr", finemap_snps$CHR),
  ranges = IRanges(start = finemap_snps$BP, end = finemap_snps$BP)
)
names(snp_ranges) <- paste0("Chr", finemap_snps$CHR, "_", finemap_snps$BP)  # Name = Chromosome_Position

# --------------------------------------------------------------------------------------------------------
# Step 6.A.2: Define GRanges object for genes (with upstream & downstream extension)
gene_ranges <- GRanges(
  seqnames = paste0("Chr", genes$chromosome_name),
  ranges = IRanges(
    start = pmax(1, genes$start_position - window), # Avoid negative values
    end = genes$end_position + window
  )
)
# Assign gene names, fallback to Ensembl ID if missing
names(gene_ranges) <- ifelse(
  genes$external_gene_name == "" | is.na(genes$external_gene_name), 
  as.character(genes$ensembl_gene_id), 
  as.character(genes$external_gene_name)
)

# --------------------------------------------------------------------------------------------------------
# Step: Find SNP-gene associations
# Step 6.A.3: Compute overlaps and extract correct names 
snp_gene_overlaps <- findOverlaps(snp_ranges, gene_ranges)

# --------------------------------------------------------------------------------------------------------
# Step 6.A.4: Extract SNP and Gene names from the matched indices
snp_gene_mapping <- data.frame(
  SNP = names(snp_ranges)[queryHits(snp_gene_overlaps)],  # Corrected extraction
  Gene = names(gene_ranges)[subjectHits(snp_gene_overlaps)]
)

# Save SNP-gene annotations
fwrite(snp_gene_mapping, file.path(output_dir, "SNP_gene_annotations.csv"), sep = "\t")
print("SNP functional annotation completed successfully!")

# --------------------------------------------------------------------------------------------------------
# Step 6.B: Credible SNPS
# Step 6.B.1: Define GRanges object for SNPs (only credible SNPs)
credible_snp_ranges <- GRanges(
  seqnames = paste0("Chr", credible_snps$CHR),
  ranges = IRanges(start = credible_snps$BP, end = credible_snps$BP)
)
names(credible_snp_ranges) <- paste0("Chr", credible_snps$CHR, "_", credible_snps$BP)  # Name = Chromosome_Position

# --------------------------------------------------------------------------------------------------------
# Step 6.B.2: Define GRanges object for genes (with upstream & downstream extension)
gene_ranges <- GRanges(
  seqnames = paste0("Chr", genes$chromosome_name),
  ranges = IRanges(
    start = pmax(1, genes$start_position - window), # Avoid negative values
    end = genes$end_position + window
  )
)
# Assign gene names, fallback to Ensembl ID if missing
names(gene_ranges) <- ifelse(
  genes$external_gene_name == "" | is.na(genes$external_gene_name), 
  as.character(genes$ensembl_gene_id), 
  as.character(genes$external_gene_name)
)

# --------------------------------------------------------------------------------------------------------
# Step: Find SNP-gene associations
# Step 6.B.3: Compute overlaps and extract correct names 
credible_snp_gene_overlaps <- findOverlaps(credible_snp_ranges, gene_ranges)

# --------------------------------------------------------------------------------------------------------
# Step 6.B.4: Extract SNP and Gene names from the matched indices
credible_snp_gene_mapping <- data.frame(
  SNP = names(credible_snp_ranges)[queryHits(credible_snp_gene_overlaps)],  # Corrected extraction
  Gene = names(gene_ranges)[subjectHits(credible_snp_gene_overlaps)]
)

# Save SNP-gene annotations
fwrite(credible_snp_gene_mapping, file.path(output_dir, "Credible_SNP_gene_annotations.csv"), sep = "\t")
print("SNP functional annotation completed successfully!")

# --------------------------------------------------------------------------------------------------------
# Step 6.B.5: Fetch gene descriptions and functional annotations from Ensembl Plants using Biomart

# Select Oryza sativa (rice) database
mart <- useMart("plants_mart", dataset = genome_db, host = "https://plants.ensembl.org")
listFilters(mart)
# Get functional annotation for the identified genes
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description", 
                 "go_id", "name_1006"),
  filters = "external_gene_name",
  values = unique(credible_snp_gene_mapping$Gene),
  mart = mart
)

print(gene_annotations)
# ensembl_gene_id external_gene_name description go_id name_1006
# 1    Os03g0407100   OSJNBa0030J19.19          NA    NA        NA
# 2    Os03g0407400                GS3          NA    NA        NA
# 3    Os03g0406900   OSJNBa0002D18.12          NA    NA        NA
# 4    Os03g0407900          OsRLCK110          NA    NA        NA

# Os03g0407100 encodes uncharacterized protein
# Os03g0407400 encodes GS3 protein that has been shown to regulate seed size in rice (Seed length and weight protein long form for short seed )
# - GS3 protein is a major regulator of grain size, specifically grain length and weight, acting as a negative regulator by controlling cell proliferation and interacting with other G-protein subunits. 
# Os03g0406900 encodes uncharacterized protein.
# Os03g0407900 encodes serine/threonine protein kinase that might have functions related to ATP binding and protein kinase activity.

# --------------------------------------------------------------------------------------------------------
# Step 7: Functional Annotation (GO & KEGG)
# ---------------------------------------------------------------------------------------------------------
# Step 7.1: Load candidate gene IDs
# These should be in RAP-DB format (e.g., Os03g0108700)
gene_ids <- genes$ensembl_gene_id

# --------------------------------------------------------------------------------------------------------
# Step 7.2: GO Enrichment Analysis (Biological Process)
go_results <- enrichGO(
  gene = gene_ids,
  OrgDb = org.Osativa.eg.db,
  keyType = "RAP",# Rice Annotation Project (not Ensembl)
  ont = "BP", # biological Process
  pAdjustMethod = "BH"
)
# enrichGO expects input gene ID to be Os06g0138000,Os11g0547000,Os05g0166900,Os02g0751800 (for rice)

# Save GO Enrichment Results
go_results_file <- file.path(output_dir, "GO_enrichment_results.csv")
fwrite(as.data.frame(go_results), go_results_file, sep = "\t")
cat("GO enrichment results saved in:", go_results_file)

# Visualize GO enrichment
dotplot(go_results, showCategory = 10, title = "GO Enrichment (BP)")

# ------------------------------------------------------------------------------------------------
# 7.3 KEGG pathway analysis 
# 7.3.1 Convert Ensembl Gene IDs to Entrez IDs using biomaRt
ensembl <- useMart("plants_mart", dataset = "osativa_eg_gene", host = "https://plants.ensembl.org")

# Retrieve Entrez Gene IDs
conversion_table <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)

# Remove NA values (genes without Entrez ID mapping)
conversion_table <- conversion_table[!is.na(conversion_table$entrezgene_id), ]

# Extract valid Entrez Gene IDs
entrez_gene_ids <- as.character(conversion_table$entrezgene_id)

# Prefix with KEGG organism code ("osa") for Oryza sativa
osa_gene_ids <- paste0("osa:", entrez_gene_ids) # Add "osa:" prefix manually

# Save KEGG conversion mapping to file
osa_conversion_df <- data.frame(Entrez_ID = entrez_gene_ids, KEGG_ID = osa_gene_ids, stringsAsFactors = FALSE)
write.table(osa_conversion_df, "kegg_manual_conversion.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# ------------------------------------------------------------------------------------------------
# Step 7.3.2: Perform KEGG pathway enrichment analysis using clusterProfiler
kegg_results <- enrichKEGG(
  gene = osa_conversion_df$Entrez_ID,
  organism = 'osa',  # KEGG code for Oryza sativa
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Display KEGG results (if any)
print(kegg_results)
# 0 enriched terms found

# ------------------------------------------------------------------------------------------------
# Step 8: Protein-Protein Interaction (PPI) Network
# ------------------------------------------------------------------------------------------------

# Load STRING database
string_db <- STRINGdb$new(version = "11.5", species = 39947, score_threshold = 400)

# ------------------------------------------------------------------------------------------------
# Step 8.1: Map gene symbols to STRING IDs
# Input: gene_vector (character vector of gene symbols)
mapped_genes <- string_db$map(data.frame(Gene = gene_vector), 
                              "Gene", 
                              removeUnmappedRows = TRUE)

# ------------------------------------------------------------------------------------------------
# Step 8.2: Retrieve PPI interactions using STRING IDs
ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)

# Save raw PPI network to file
if (!is.null(ppi_network)) {
  fwrite(ppi_network, file.path(output_dir, "PPI_network.csv"), sep = "\t")
  cat("PPI network analysis completed and saved!")
} else {
  cat("No significant protein-protein interactions found.")
}

# ------------------------------------------------------------------------------------------------
# Step 8.3: Clean & FIlter PPI Network
# - Remove duplicate interactions
# - Filter by confidence score (> 700 indicates high confidence)
ppi_network <- unique(ppi_network)
ppi_network <- ppi_network[ppi_network$combined_score > 700, ]

# ------------------------------------------------------------------------------------------------
# Step 8.4: Network Visualization using igraph and ggraph
library(igraph)
library(ggraph)
library(tidygraph)

# Convert interaction table into an igraph object
ppi_graph <- graph_from_data_frame(ppi_network, directed = FALSE)

# Plot the PPI Network
ggraph(ppi_graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = combined_score), show.legend = FALSE) +
  geom_node_point(size = 5, color = "blue") +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()

# ------------------------------------------------------------------------------------------------
# Step 8.5: Topological Network Analysis
# Analyze node properties to understand key genes in the PPI network

# A. Degree centrality (number of connections)
degree_centrality <- degree(ppi_graph)
top_degrees <- sort(degree_centrality, decreasing = TRUE)[1:10]
print(top_degrees)

# B. Betweenness centrality (nodes that act as bridges)
betweenness_centrality <- betweenness(ppi_graph)

# C. Clustering Coefficient (how tightly clustered each node's neighbors are)
clustering_coeff <- transitivity(ppi_graph, type = "local")

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------