# Since these steps are repeated for all the different analysis scripts, the decision was made to separate 
# them into a different scrpt that can be loaded useing source(filename) at any given point. 
# This script is used to : 
# 1. Construct a DESeqDataSet object with :
#    - raw counts, 
#    - sample metadata, and 
#    - a design formula to model differences based on the group variable.
# 2. Enables inspection of the DESeq2 dataset for differential expression analysis.


# Important to set the working directory as the one holfing the reformatted_counts.txt. 
# Otherwise the code will not work 
# 1. install.packages("BiocManager")
# 2. BiocManager::install("DESeq2")
# 3. install.packages("ggplot2")
# 4. install.packages("pheatmap")

# Importing the necessary libabries, both of them had to be downloaded first using the following two commands 
# loading DeSeq2
library(DESeq2)
# loading ggplot2
library(ggplot2)
# loading pheatmap 
library(pheatmap)
package_version <- packageVersion("pheatmap")
cat("clusterProfiler version:", as.character(package_version), "\n")


# Read the formatted featureCounts file
# row.names = 1 tells R to treat the first column (gene names) as the row names
# meaning they are no longer part of the main data matrix 
counts <- read.table("reformatted_counts.txt", header = TRUE, row.names = 1)


# adds the appropriate column names to the data matrix meaning that now both 
# rows and columns can be appropriately identified 
colnames(counts) <- c(
  "SRR7821951", "SRR7821969", "SRR7821920", "SRR7821938", 
  "SRR7821922", "SRR7821939", "SRR7821918", "SRR7821968", 
  "SRR7821953", "SRR7821949", "SRR7821970", "SRR7821950", 
  "SRR7821952", "SRR7821919", "SRR7821937", "SRR7821921"
)

# creates a data frame (sampleInfo) providing metadata about the samples 
# rename the sample 
sampleInfo <- data.frame(
  row.names = colnames(counts),
  # name the samples according to the information provided by the README file in 
  # the directory with the original 
  group = c(
    "Blood_Case", "Blood_Control", "Lung_Case", "Lung_Control",
    "Lung_Case", "Lung_Control", "Lung_Case", "Blood_Control",
    "Blood_Case", "Blood_Case", "Blood_Control", "Blood_Case",
    "Blood_Case", "Lung_Case", "Lung_Control", "Lung_Case"
  )
)
# column 1 : Sample Names ("SRR7821951",...)
# column 2 : Corresponding groups : (Blood_Case,...) 


# Create the DESeqDataSet object to store the necessary information 
dds <- DESeqDataSetFromMatrix(
  # raw count data matrix, rows : genes
  countData = counts,
  # metadata : samples corresponding groups 
  colData = sampleInfo,
  # tilde : specifies model formula
  # group is the factor used to model differences in gene expression between conditions
  design = ~ group
)

# to inspect the DESeqDataSet object, make sure the values are in a realistic range 
dds
# To make sure that it has creted the object correctly, we ask the following question : 
# 1. DESeqDataSet object has the expected dimensions (genes × samples)? 

# Run DESeq normalisation
# - Adjusts raw counts for sequencing depth and other biases
#   - Sequencing depth : total number of RNA-sequencing reads obtained for a sample
#   - If Sample A has 10 million reads, and a gene has 1,000 counts, it means the gene accounts for 0.01% of total reads

# Dispersion Estimation: Estimates variability in gene expression across replicates.
# Model Fitting: Fits a negative binomial model to the count data for differential expression analysis.
# Result: dds now contains normalized counts and fitted models, ready for downstream analysis.
dds <- DESeq(dds)

# Avaliable values after this step : 
#     - normalised counts : counts(dds, normalized = TRUE)
#     - Dispersion Estimates (per-gene dispersion (α)) : dispersions(dds)
#     - Size factors : sizeFactors(dds) (used to normalise each count)
#     - model coefficients : coef(dds) 
#     - Fitted Values: assays(dds)$mu (expected count values for each gene)
#     - A result table : results_table <- results(dds) including : 
#           - log2 fold change (LFC)
#           - Standard error (SE) for LFC
#           - Wald statistic
#           - P-value (for LFC) 
#           - Adjusted p-value (padj) (corrected with BH correction)
















