# GWAS_rice: Genomic Prediction and Fine-Mapping in Rice

## Overview

This repository presents a complete genome-wide association study (GWAS) and trait prediction pipeline applied to rice (*Oryza sativa*) to identify and interpret genomic variants associated with grain weight. It integrates:

- Phenotypic and genotypic preprocessing  
- GWAS using the GAPIT package  
- Bayesian fine-mapping and candidate gene identification  
- Functional enrichment (GO, KEGG, PPI)  
- Machine learning models for genomic prediction

## Dataset Information

- **Phenotype and genotype source:**  
  - https://www.ebi.ac.uk/ena/browser/view/PRJEB6064  
  - https://iagr.genomics.cn/CropGS/#/Datasets  
- **Genotype format:** PLINK binary  
- **Traits measured:** Agronomic traits from multi-location rice field trials  
- **Target trait for prediction:** Grain weight

---

## Workflow Overview

The project is organized into eight major steps:

### Step 1: Phenotypic Data Processing (`Phenotypic_data_processing.R`)

- Computes Best Linear Unbiased Predictions (BLUPs) to account for location effects  
- Performs quality control, normality testing, and Inverse Normal Transformation (INT)  
- Outputs a cleaned phenotype file for GWAS and prediction models

**Key outputs:**

- `blup_results/computed_BLUP_from_SY_HZ.csv`  
- `GAPIT_data/pheno_BLUP_INT_cleaned.txt`  
- Visualizations: Q-Q plots, histograms, correlation heatmaps

---

### Step 2: Genotypic Data Processing (Bash Scripts)

Scripts:  
- `genotype_qc_pipeline.sh`  
- `fix_sample_id_in_vcf.sh`  
- `convert_vcf_to_hapmap.sh`

Processing steps:  
- Quality control using PLINK (MAF, missingness)  
- VCF header fixing using `bcftools`  
- VCF-to-HapMap conversion using TASSEL

**Final output:**  
- `gwas_ready.hmp.txt` (ready for GAPIT)

---

### Step 3: GWAS Analysis using GAPIT (`GAPIT_analysis.R`)

- GWAS performed using the FarmCPU model  
- Controls for population structure using 3 principal components  
- Identifies SNPs significantly associated with grain weight

**Outputs:**

- `.GWAS.Results.csv` file with p-values  
- Manhattan and Q-Q plots  
- PCA structure matrix

---

### Step 4: Post-GWAS Fine-Mapping and Candidate Gene Prioritization

Scripts:  
- `pGWAS_01_extract_alleles.R`  
- `pGWAS_02_clumping.sh`  
- `pGWAS_03_finemap.sh`  
- `pGWAS_03b_finemap.R`  
- `pGWAS_04_candidate_genes_final.R`

Steps:  
- PLINK-based clumping to select lead SNPs  
- Bayesian fine-mapping using SuSiE to identify credible SNPs (PIP > 0.95)  
- Gene annotation with biomaRt  
- GO and KEGG enrichment using clusterProfiler  
- Protein-protein interaction (PPI) network with STRINGdb and igraph

**Key outputs:**

- `credibleSNPs_candidate_genes.csv`  
- `GO_enrichment_results.csv`  
- `PPI_network.csv`

---

### Step 5: Genomic Prediction — G-BLUP and Bayesian Models

- Computed the genomic relationship matrix (G-matrix) using GWAS-significant SNPs  
- Performed G-BLUP using:  
  - Ridge regression  
  - Kernel Ridge regression  
- Performed Bayesian regression using cmdstanpy:  
  - Bayesian Ridge, BayesA, BayesB, Bayesian LASSO

| Model           | RMSE   | R²     |
|-----------------|--------|--------|
| Bayesian Ridge  | 0.819  | 0.3342 |
| BayesA          | 0.9791 | 0.0484 |
| BayesB          | 0.8742 | 0.2414 |
| Bayesian LASSO  | 0.8351 | 0.3076 |

Bayesian Ridge performed best, while BayesA showed poor generalization — possibly due to overfitting.

---

### Step 6: Machine Learning Models — Random Forest and XGBoost

- Used GWAS-significant SNPs as features  
- Hyperparameters tuned using RandomizedSearchCV (5-fold CV)

| Model         | RMSE   | R²     |
|---------------|--------|--------|
| Random Forest | 0.7757 | 0.4027 |
| XGBoost       | 0.7596 | 0.4272 |

XGBoost achieved the highest predictive accuracy overall.

---

### Step 7: SNP-SNP Interaction and Feature Selection

**Interaction modeling:**

- Selected top 50 SNPs (based on feature importance from RF + XGBoost)  
- Generated pairwise SNP interactions using sklearn PolynomialFeatures  

| Model              | RMSE   | R²     |
|--------------------|--------|--------|
| RF + Interactions  | 0.7835 | 0.3906 |
| XGB + Interactions | 0.7801 | 0.3960 |

**Feature selection with top 100 SNPs:**

| Model         | RMSE   | R²     |
|---------------|--------|--------|
| RF (Top SNPs) | 0.7739 | 0.4054 |
| XGB (Top SNPs)| 0.7628 | 0.4224 |

Feature selection improved interpretability without sacrificing performance.

---

### Step 8: SNP Encoding and SVM Modeling

**One-hot encoding of SNPs:**

| Model         | RMSE   | R²     |
|---------------|--------|--------|
| RF (One-Hot)  | 0.7780 | 0.3991 |
| XGB (One-Hot) | 0.7660 | 0.4176 |

**Support Vector Regression (SVM):**

| Model         | RMSE   | R²     |
|---------------|--------|--------|
| SVM (RBF)     | 0.7892 | 0.3817 |

SVM performed well but was slower to train than tree-based models.

---

## Repository Structure

```
GWAS_rice/  
│── data/                                  # Cleaned phenotype and genotype data  
│   │── pheno_BLUP_INT_cleaned.csv         # INT-transformed phenotype file  
│   │── gwas_ready.hmp.txt                 # HapMap-formatted genotype file for GAPIT  
│  
│── scripts/                               # All analysis scripts grouped by language  
│   │  
│   ├── R/                                 # R scripts for phenotype processing, GWAS, fine-mapping  
│   │   │── Phenotypic_data_processing.R         # BLUP estimation, INT, and trait QC  
│   │   │── GAPIT_analysis.R                     # GWAS using FarmCPU in GAPIT  
│   │   │── pGWAS_01_extract_alleles.R           # Effect size calculation and allele formatting  
│   │   │── pGWAS_03b_finemap.R                  # Bayesian fine-mapping using SuSiE  
│   │   │── pGWAS_04_candidate_genes_final.R     # Candidate gene annotation and enrichment analysis  
│   │  
│   ├── bash/                              # Bash scripts for genotype QC and formatting  
│   │   │── 1_genotype_qc_pipeline.sh              # PLINK-based filtering and sample matching  
│   │   │── 2_fix_sample_id_in_vcf.sh              # bcftools-based VCF header correction  
│   │   │── 3_convert_vcf_to_hapmap.sh             # TASSEL-based VCF to HapMap conversion  
│   │   │── pGWAS_02_clumping.sh                 # PLINK-based LD clumping of SNPs  
│   │   │── pGWAS_03_finemap.sh                  # Prepares LD matrix and FINEMAP input files  
│   │  
│   ├── python/                            # Python-based machine learning models and predictions  
│   │   │── GWAS_ML_Single_Trait_Prediction.ipynb  # Genomic prediction using ML and Bayesian models  
│  
│── models/                                # Trained models and Stan model files  
│   │── stan/                              # Stan models for BayesA, BayesB, Ridge, LASSO  
│   │── trained_ml_models.pkl              # Trained Random Forest / XGBoost models  
│  
│── results/                                # Organized outputs from the entire pipeline 
│   │── blup_results/                        # BLUP estimation, normality, INT plots and CSVs
│   │── GAPIT_data/                         # GWAS summary statistics, plots, PCA 
│   │── genotype_qc_analysis/                # Filtered genotypes, VCFs, HapMap file                  
│   │── GWAS_GAPIT/                         # GAPIT-generated files (.csv, .png, .pdf) 
│   │── Machine_learning/                   # ML & Bayesian model performance scores, predictions 
│   │── pGWAS/                              # Fine-mapping outputs, credible SNPs, annotations  
│  
│── README.md                              # Project documentation  

```
---

## Running the Workflow

You can run each part in sequence:

1. **Phenotype preprocessing:**  
   `Rscript scripts/Phenotypic_data_processing.R`

2. **Genotype QC and conversion:**  
   Run shell scripts from `bash/` in order

3. **GWAS in R:**  
   `Rscript scripts/GAPIT_analysis.R`

4. **Fine-mapping & Annotation:**  
   Run `pGWAS_*` scripts as per the order listed

5. **Genomic prediction models:**  
   Use Python scripts or notebooks for training and evaluation

---

## Author

**Abhishek Shrestha, Ph.D**  
Email: abhishek.shrestha39@gmail.com  
Location: Mesa, Arizona