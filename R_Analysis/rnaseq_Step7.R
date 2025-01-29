#------------------------ Libraries and Setup ----------------------------------
  
# Load necessary libraries
  
# To perform GO enrichment analysis
library(clusterProfiler)
# Annotation database for mouse
library(org.Mm.eg.db)
# For additional visualisation and analysis options
library(DOSE)
# For visualisations
library(ggplot2)
# to place several plots next to each other 
library(gridExtra)

# ------------------------ Data Preparation ------------------------------------
  
# Load data from the previous file (Step 6)
source("rnaseq_step6.R")

# Define diff. Expressed Genes

# Genes extracted based on adjusted p-value < 0.05 for blood
# and lung data sets separately

diff_expr_genes_blood <- rownames(res_blood[res_blood$padj < 0.05 & !is.na(res_blood$padj), ])
diff_expr_genes_lung <- rownames(res_lung[res_lung$padj < 0.05 & !is.na(res_lung$padj), ])

# Combine DEGs from both comparisons
diff_expr_genes_combined <- unique(c(diff_expr_genes_blood, diff_expr_genes_lung))

# Define the universe of all measured genes
# Contains the Ensembl IDs of all genes
all_genes <- rownames(counts)

# Ensure all genes have been loaded
print(length(all_genes) == nrow(counts))

#------------------------ GO Enrichment Analysis: Blood ------------------------
  
# Perform GO enrichment analysis for blood
  
  enrich_res_blood <- enrichGO(
    gene = diff_expr_genes_blood,          # List of different genes 
    universe = all_genes,                  # Universe: all relevant genes
    OrgDb = org.Mm.eg.db,                  # Annotation database
    ont = "ALL",                           # GO subontology: BP, MF, CC, or ALL
    keyType = "ENSEMBL",                   # Gene ID type: ENSEMBL IDs
    pvalueCutoff = 0.05,                   # Adjusted p-value cutoff
    qvalueCutoff = 0.2,                    # False discovery rate cutoff
    readable = TRUE                        # Convert gene IDs to gene names
  )

# Display top results in the console for blood
cat("Top GO terms for blood by category:\n")
print(head(as.data.frame(enrich_res_blood)))

#------------------------ GO Enrichment Analysis: Lung -----------------
  
  # Perform GO enrichment analysis for lung
  
  enrich_res_lung <- enrichGO(
    gene = diff_expr_genes_lung,          # List of DEGs
    universe = all_genes,                 # Universe: all genes measured
    OrgDb = org.Mm.eg.db,                 # Annotation database
    ont = "ALL",                          # GO subontology: BP, MF, CC, or ALL
    keyType = "ENSEMBL",                  # Gene ID type: ENSEMBL IDs
    pvalueCutoff = 0.05,                  # Adjusted p-value cutoff
    qvalueCutoff = 0.2,                   # False discovery rate cutoff
    readable = TRUE                       # Convert gene IDs to gene names
  )

# Display top results in the console for lung
cat("Top GO terms for lung by category:\n")
print(head(as.data.frame(enrich_res_lung)))


#------------------------ GO Enrichment Analysis: Combined ---------------------
  
  # Perform GO enrichment analysis for combined DEGs
  
  enrich_res_combined <- enrichGO(
    gene = diff_expr_genes_combined,      # List of combined DEGs
    universe = all_genes,                 # Universe: all genes measured
    OrgDb = org.Mm.eg.db,                 # Annotation database
    ont = "ALL",                          # GO subontology: BP, MF, CC, or ALL
    keyType = "ENSEMBL",                  # Gene ID type: ENSEMBL IDs
    pvalueCutoff = 0.05,                  # Adjusted p-value cutoff
    qvalueCutoff = 0.2,                   # False discovery rate cutoff
    readable = TRUE                       # Convert gene IDs to gene names
  )

# Display top results in the console for combined DEGs
cat("Top GO terms for combined DEGs by category:\n")
print(head(as.data.frame(enrich_res_combined)))


#------------------------ Combine Plots ----------------------------------------
  
# Arrange dot plots for blood, lung, and combined next to each other

# create different options to pick from for the report later on 

# Generate dot plots for each dataset
dotplot_blood <- dotplot(enrich_res_blood, showCategory = 10, title = "GO Term Enrichment (Blood)") + 
  theme(plot.title = element_text(hjust = 0.5))

dotplot_lung <- dotplot(enrich_res_lung, showCategory = 10, title = "GO Term Enrichment (Lung)") + 
  theme(plot.title = element_text(hjust = 0.5))

dotplot_combined <- dotplot(enrich_res_combined, showCategory = 10, title = "GO Term Enrichment (Combined)") + 
  theme(plot.title = element_text(hjust = 0.5))

# Arrange plots side by side
grid.arrange(
  dotplot_blood,
  dotplot_lung,
  dotplot_combined,
  ncol = 3
)

# Arrange plots in two columns
grid.arrange(
  dotplot_blood,
  dotplot_lung,
  ncol = 2
)

# ------------------- Analysis on overlapping genes ----------------------------


# Identify overlapping DEGs between blood and lung
overlapping_genes <- intersect(diff_expr_genes_blood, diff_expr_genes_lung)

# Check how many genes overlap
length(overlapping_genes)

# Print overlapping genes
print(overlapping_genes)

# Perform GO enrichment analysis for overlapping genes
enrich_res_overlap <- enrichGO(
  gene = overlapping_genes,      # Overlapping DEGs
  universe = all_genes,          # Background: all measured genes
  OrgDb = org.Mm.eg.db,          # Annotation database
  ont = "ALL",                   # Use all GO subontologies
  keyType = "ENSEMBL",           # Gene ID type: ENSEMBL
  pvalueCutoff = 0.05,           # Adjusted p-value cutoff
  qvalueCutoff = 0.2,            # FDR cutoff
  readable = TRUE                # Convert IDs to gene names
)

library(pheatmap)


# Extract expression data for overlapping genes
overlap_expr_data <- counts[overlapping_genes, ]

# Calculate variance for each gene
gene_variances <- apply(overlap_expr_data, 1, var)

# Select top 30 most variable genes
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:30]

# Filter the expression data for top 100 genes
top_gene_expr_data <- overlap_expr_data[top_genes, ]

# Normalise counts 
log_counts <- log2(top_gene_expr_data + 1)

# Generate heatmap for top 100 genes
pheatmap(
  log_counts,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,  # Show row names for top genes
  fontsize_row = 6       # Adjust font size for row names
)




