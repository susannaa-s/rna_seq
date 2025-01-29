# STEP 6: Differential Expression Analysis
# This script performs differential expression analysis to identify significant genes
# and visualize the results using volcano plots for two comparisons: Blood and Lung.

# Load required librariesy
library(ggrepel)  # For adding non-overlapping text labels to plots
library(pheatmap)  # For heatmaps (not used here but often in differential analysis workflows)
library(dplyr)  # For data manipulation
library(patchwork)  # For combining plots
library(tibble)  # For adding row names as a column

# Load the custom function to create the DESeqDataSet object
source("create_DESeqDataSet.R")

# --------------------------- CASE vs. CONTROL --------------------------------------------

# Extract differential expression results
# DESeq2 produces DESeqResults objects that contain information about log2 fold changes,
# p-values, adjusted p-values (padj), and other metrics.
# The contrast argument specifies the comparison of interest:
#     - Blood: Case vs. Control
#     - Lung: Case vs. Control
res_blood <- results(dds, contrast = c("group", "Blood_Case", "Blood_Control"))
res_lung <- results(dds, contrast = c("group", "Lung_Case", "Lung_Control"))

# Identify the top 5 most significant genes based on adjusted p-values (padj)
# Sorting by padj (ascending) ensures the most statistically significant genes are selected.
top_blood <- res_blood[order(res_blood$padj), ][1:5, ]  # Top 5 genes for Blood
top_lung <- res_lung[order(res_lung$padj), ][1:5, ]  # Top 5 genes for Lung

# Extract gene names for the top genes
top_blood_genes <- rownames(top_blood)  # Gene names for Blood
top_lung_genes <- rownames(top_lung)  # Gene names for Lung

# Create data frames with annotations for top genes
# Add gene names as a column and assign short labels (e.g., G1, G2) for plot annotations.
top_blood_genes_df <- as.data.frame(res_blood) %>%
  rownames_to_column(var = "gene") %>%
  filter(gene %in% top_blood_genes) %>%
  mutate(short_label = paste0("G", seq_len(nrow(.))))  # Short labels for top genes

top_lung_genes_df <- as.data.frame(res_lung) %>%
  rownames_to_column(var = "gene") %>%
  filter(gene %in% top_lung_genes) %>%
  mutate(short_label = paste0("G", seq_len(nrow(.))))  # Short labels for top genes

# Print the data frames for verification
print(top_blood_genes_df)
print(top_lung_genes_df)

# Set the significance threshold
alpha <- 0.05  # Genes with adjusted p-values < alpha are considered significant


# Process Blood data for plotting
res_blood_df <- as.data.frame(res_blood) %>%
  filter(!is.na(padj))  # Exclude genes with NA adjusted p-values

# Handle cases where padj is 0 by replacing with a small value (to avoid -Inf in log10)
min_padj_blood <- min(res_blood_df$padj[res_blood_df$padj > 0], na.rm = TRUE)  # Smallest non-zero padj
res_blood_df <- res_blood_df %>%
  mutate(
    plot_padj = -log10(ifelse(padj == 0, min_padj_blood / 10, padj)),  # Adjusted p-value for plotting
    # Categorize genes based on statistical and biological significance
    point_category = case_when(
      padj < alpha & abs(log2FoldChange) >= 1 ~ "High significance",
      padj < alpha & abs(log2FoldChange) < 1 ~ "Low biological impact",
      TRUE ~ "Not significant"
    )
  )
max_y_blood <- max(res_blood_df$plot_padj, na.rm = TRUE)  # Maximum y-axis value for Blood plot

# Volcano plot for Blood
volcano_blood <- ggplot(data = res_blood_df, aes(x = log2FoldChange, y = plot_padj)) +
  geom_point(aes(color = point_category), alpha = 0.6) +  # Plot points with colors by significance
  scale_color_manual(values = c(
    "High significance" = "deepskyblue4",  # Significant genes with large effect size in deep blue
    "Low biological impact" = "lightblue",  # Significant genes with small effect size in light blue
    "Not significant" = "grey"  # Non-significant genes in grey
  )) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkslategrey") +  # Threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +  # Effect size thresholds
  scale_y_continuous(limits = c(0, max_y_blood + 10), expand = c(0, 0)) +  # Adjust y-axis range
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Blood : Case vs. Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Category"
  ) +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

print(volcano_blood)

# Process Lung data for plotting
res_lung_df <- as.data.frame(res_lung) %>%
  filter(!is.na(padj))  # Exclude genes with NA adjusted p-values

# Handle cases where padj is 0 by replacing with a small value (to avoid -Inf in log10)
min_padj_lung <- min(res_lung_df$padj[res_lung_df$padj > 0], na.rm = TRUE)  # Smallest non-zero padj
res_lung_df <- res_lung_df %>%
  mutate(
    plot_padj = -log10(ifelse(padj == 0, min_padj_lung / 10, padj)),  # Adjusted p-value for plotting
    # Categorize genes based on statistical and biological significance
    point_category = case_when(
      padj < alpha & abs(log2FoldChange) >= 1 ~ "High significance",
      padj < alpha & abs(log2FoldChange) < 1 ~ "Low biological impact",
      TRUE ~ "Not significant"
    )
  )
max_y_lung <- max(res_lung_df$plot_padj, na.rm = TRUE)  # Maximum y-axis value for Lung plot

# Volcano plot for Lung
volcano_lung <- ggplot(data = res_lung_df, aes(x = log2FoldChange, y = plot_padj)) +
  geom_point(aes(color = point_category), alpha = 0.6) +
  scale_color_manual(values = c(
    "High significance" = "darkgreen",  # Significant genes with large effect size in dark green
    "Low biological impact" = "lightgreen",  # Significant genes with small effect size in light green
    "Not significant" = "grey"  # Non-significant genes in grey
  )) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkslategrey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +
  scale_y_continuous(limits = c(0, max_y_lung + 10), expand = c(0, 0)) +
  theme_minimal() +
  labs(
    title = "Lung : Case vs. Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Category"
  ) +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

print(volcano_lung)

# Combine and display the plots
combined_plot <- volcano_blood + volcano_lung
print(combined_plot)



# ------------------------------ Blood vs Lung ----------------------------------

# This section is not mentionned or discussed in the report. Served as additional information to
# to look at the data from a different perspective as well.  
# Set the significance threshold
alpha <- 0.05  # Genes with adjusted p-values < alpha are considered significant

# Process Lung_Case vs Blood_Case data for plotting
res_bl_case_df <- as.data.frame(results(dds, contrast = c("group", "Lung_Case", "Blood_Case"))) %>%
  filter(!is.na(padj))  # Exclude genes with NA adjusted p-values

# Handle cases where padj is 0 by replacing with a small value (to avoid -Inf in log10)
min_padj <- min(res_bl_case_df$padj[res_bl_case_df$padj > 0], na.rm = TRUE)  # Smallest non-zero padj
res_bl_case_df <- res_bl_case_df %>%
  mutate(
    plot_padj = -log10(ifelse(padj == 0, min_padj / 10, padj)),  # Adjusted p-value for plotting
    # Categorize genes based on statistical and biological significance
    point_category = case_when(
      padj < alpha & abs(log2FoldChange) >= 1 ~ "High significance",
      padj < alpha & abs(log2FoldChange) < 1 ~ "Low biological impact",
      TRUE ~ "Not significant"
    )
  )
max_y_bl_case <- max(res_bl_case_df$plot_padj, na.rm = TRUE)  # Maximum y-axis value for this plot

# Volcano plot for Lung_Case vs Blood_Case
volcano_bl_case <- ggplot(data = res_bl_case_df, aes(x = log2FoldChange, y = plot_padj)) +
  geom_point(aes(color = point_category), alpha = 0.6) +  # Plot points with colors by significance
  scale_color_manual(values = c(
    "High significance" = "darkblue",  # Significant genes with large effect size in dark blue
    "Low biological impact" = "lightblue",  # Significant genes with small effect size in light blue
    "Not significant" = "grey"  # Non-significant genes in grey
  )) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkslategrey") +  # Threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +  # Effect size thresholds
  scale_y_continuous(limits = c(0, max_y_bl_case + 10), expand = c(0, 0)) +  # Adjust y-axis range
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Lung vs. Blood : Case",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Category"
  ) +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

print(volcano_bl_case)

# Process Lung_Control vs Blood_Control data for plotting
res_bl_control_df <- as.data.frame(results(dds, contrast = c("group", "Lung_Control", "Blood_Control"))) %>%
  filter(!is.na(padj))  # Exclude genes with NA adjusted p-values

# Handle cases where padj is 0 by replacing with a small value (to avoid -Inf in log10)
min_padj <- min(res_bl_control_df$padj[res_bl_control_df$padj > 0], na.rm = TRUE)  # Smallest non-zero padj
res_bl_control_df <- res_bl_control_df %>%
  mutate(
    plot_padj = -log10(ifelse(padj == 0, min_padj / 10, padj)),  # Adjusted p-value for plotting
    # Categorize genes based on statistical and biological significance
    point_category = case_when(
      padj < alpha & abs(log2FoldChange) >= 1 ~ "High significance",
      padj < alpha & abs(log2FoldChange) < 1 ~ "Low biological impact",
      TRUE ~ "Not significant"
    )
  )
max_y_bl_control <- max(res_bl_control_df$plot_padj, na.rm = TRUE)  # Maximum y-axis value for this plot

# Volcano plot for Lung_Control vs Blood_Control
volcano_bl_control <- ggplot(data = res_bl_control_df, aes(x = log2FoldChange, y = plot_padj)) +
  geom_point(aes(color = point_category), alpha = 0.6) +
  scale_color_manual(values = c(
    "High significance" = "darkgreen",  # Significant genes with large effect size in dark green
    "Low biological impact" = "lightgreen",  # Significant genes with small effect size in light green
    "Not significant" = "grey"  # Non-significant genes in grey
  )) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkslategrey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +
  scale_y_continuous(limits = c(0, max_y_bl_control + 10), expand = c(0, 0)) +
  theme_minimal() +
  labs(
    title = "Lung vs. Blood : Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Category"
  ) +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

print(volcano_bl_control)

# Combine and display the plots
combined_cross_tissue_plot <- volcano_bl_case + volcano_bl_control
print(combined_cross_tissue_plot)


# -------------------------  Questions -----------------------------------------

# Set the adjusted p-value threshold for significance
alpha <- 0.05  # Genes with adjusted p-values < 0.05 are considered significant

# Q1: How many genes are differentially expressed (DE) in the pairwise comparison?
cat("\nNumber of differentially expressed genes (DE):\n")
cat("For Blood:\n")
cat(sum(res_blood$padj < alpha, na.rm = TRUE), "genes\n")
cat("For Lung:\n")
cat(sum(res_lung$padj < alpha, na.rm = TRUE), "genes\n")
cat("For Case (Lung_Case vs Blood_Case):\n")
cat(sum(res_bl_case_df$padj < alpha, na.rm = TRUE), "genes\n")
cat("For Control (Lung_Control vs Blood_Control):\n")
cat(sum(res_bl_control_df$padj < alpha, na.rm = TRUE), "genes\n")

# Q2: How many DE genes are up-regulated vs down-regulated?
cat("\nUpregulated vs Downregulated DE genes:\n")
cat("For Blood:\n")
cat(sum(res_blood$padj < alpha & res_blood$log2FoldChange > 0, na.rm = TRUE), "upregulated genes\n")
cat(sum(res_blood$padj < alpha & res_blood$log2FoldChange < 0, na.rm = TRUE), "downregulated genes\n")
cat("For Lung:\n")
cat(sum(res_lung$padj < alpha & res_lung$log2FoldChange > 0, na.rm = TRUE), "upregulated genes\n")
cat(sum(res_lung$padj < alpha & res_lung$log2FoldChange < 0, na.rm = TRUE), "downregulated genes\n")
cat("For Case (Lung_Case vs Blood_Case):\n")
cat(sum(res_bl_case_df$padj < alpha & res_bl_case_df$log2FoldChange > 0, na.rm = TRUE), "upregulated genes\n")
cat(sum(res_bl_case_df$padj < alpha & res_bl_case_df$log2FoldChange < 0, na.rm = TRUE), "downregulated genes\n")
cat("For Control (Lung_Control vs Blood_Control):\n")
cat(sum(res_bl_control_df$padj < alpha & res_bl_control_df$log2FoldChange > 0, na.rm = TRUE), "upregulated genes\n")
cat(sum(res_bl_control_df$padj < alpha & res_bl_control_df$log2FoldChange < 0, na.rm = TRUE), "downregulated genes\n")


# Q3 : Based on the original publication, select 2-3 genes that are of particular interest and investigate their expression level.

# ----------------------------- Setup and Preprocessing ---------------------------------

# Retrieve normalised counts from DESeqDataSet
norm_counts <- counts(dds, normalized = TRUE)

# Define Sample Metadata
sampleInfo <- data.frame(
  row.names = colnames(norm_counts),
  Sample = colnames(norm_counts),
  Condition = c(
    "Blood_Case", "Blood_Control", "Lung_Case", "Lung_Control",
    "Lung_Case", "Lung_Control", "Lung_Case", "Blood_Control",
    "Blood_Case", "Blood_Case", "Blood_Control", "Blood_Case",
    "Blood_Case", "Lung_Case", "Lung_Control", "Lung_Case"
  )
)

# Define Genes of Interest
selected_genes <- c(
  "ENSMUSG00000000386",  # Mx1
  "ENSMUSG00000028270",  # Gbp2
  "ENSMUSG00000074896",  # Ifit3 
  "ENSMUSG00000018899"   # Irf1
)

# Convert res_blood and res_lung to data frames
res_blood_df <- as.data.frame(res_blood) %>%
  rownames_to_column(var = "Gene") %>%
  mutate(Source = "Blood")

res_lung_df <- as.data.frame(res_lung) %>%
  rownames_to_column(var = "Gene") %>%
  mutate(Source = "Lung")

# Combine normalised counts with metadata
norm_counts_long <- as.data.frame(norm_counts) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(
    cols = starts_with("SRR"),  # Match sample names in norm_counts
    names_to = "Sample",        # Each column becomes a row
    values_to = "Expression"    # Gene expression values
  ) %>%
  left_join(sampleInfo, by = "Sample")  # Add condition-specific metadata

# Merge gene-level data with condition-specific metadata
res_blood_df <- res_blood_df %>%
  left_join(norm_counts_long, by = "Gene") %>%
  mutate(Source = paste(Source, Condition, sep = "_")) %>%  # Combine Source and Condition
  select(-Condition)  # Remove Condition column

res_lung_df <- res_lung_df %>%
  left_join(norm_counts_long, by = "Gene") %>%
  mutate(Source = paste(Source, Condition, sep = "_")) %>%  # Combine Source and Condition
  select(-Condition)  # Remove Condition column

# Combine blood and lung results for global analysis
res_global <- bind_rows(res_blood_df, res_lung_df)

# ------------------------- Adjusted P-Value Based Analysis -------------------------

# Deduplicate to ensure each gene is ranked only once globally
# For each gene, keep the minimum adjusted p-value (padj) and the corresponding log2 fold change
res_global_unique_pval <- res_global %>%
  group_by(Gene) %>%
  reframe(
    baseMean = mean(baseMean, na.rm = TRUE),               # Average baseMean across all samples
    log2FoldChange = if (all(is.na(padj))) NA else log2FoldChange[which.min(padj)], # Log2 FC at the minimum padj
    padj = if (all(is.na(padj))) NA else min(padj, na.rm = TRUE)  # Minimum adjusted p-value for each gene
  ) %>%
  filter(!is.na(padj)) %>%                                 # Remove genes with no valid padj
  arrange(padj) %>%                                        # Sort by adjusted p-value
  mutate(GlobalRank = row_number())                        # Assign global rank based on padj

# Filter preselected genes and include their global ranks
selected_genes_global_pval <- res_global_unique_pval %>%
  filter(Gene %in% selected_genes) %>%
  arrange(GlobalRank) # Sort by global rank for clarity

# Print the selected genes with their global ranks for reporting
print(selected_genes_global_pval)

# Save the summarised results for preselected genes (optional)
# write.csv(selected_genes_global_pval, "selected_genes_global_pval.csv", row.names = FALSE)

# Extract and print the top 20 genes globally by adjusted p-value for context
top20_global_pval <- res_global_unique_pval %>%
  slice_head(n = 20)
print(top20_global_pval)

# ------------------------- Log2 Fold Change Based Analysis -------------------------

# Deduplicate to ensure each gene is ranked only once globally
# Handle NA values properly and remove them before using which.max()
res_global_unique_fc <- res_global %>%
  group_by(Gene) %>%
  reframe(
    baseMean = mean(baseMean, na.rm = TRUE),                                # Average baseMean across all samples
    log2FoldChange = if (all(is.na(log2FoldChange))) NA else log2FoldChange[which.max(abs(log2FoldChange[!is.na(log2FoldChange)]))], # Max absolute log2 FC
    padj = if (all(is.na(padj))) NA else min(padj, na.rm = TRUE)            # Minimum adjusted p-value
  ) %>%
  filter(!is.na(log2FoldChange)) %>%                                        # Remove entries with no valid log2FoldChange
  arrange(desc(abs(log2FoldChange))) %>%                                    # Sort by absolute log2 fold change
  mutate(GlobalRank = row_number())                                         # Assign global rank


# Filter preselected genes and include their global ranks
selected_genes_global_fc <- res_global_unique_fc %>%
  filter(Gene %in% selected_genes) %>%
  arrange(GlobalRank) # Sort by global rank for clarity

# Print the selected genes with their global ranks for reporting
print(selected_genes_global_fc)

# Save the summarised results for preselected genes (optional)
# write.csv(selected_genes_global_fc, "selected_genes_global_fc.csv", row.names = FALSE)

# Extract and print the top 20 genes globally by absolute log2 fold change for context
top20_global_fc <- res_global_unique_fc %>%
  slice_head(n = 20)
print(top20_global_fc)

# ------------------------- Dataset Overview and Rankings Context -------------------------

# Calculate the total number of unique genes in the dataset
total_genes <- nrow(res_global %>% distinct(Gene))
cat("Total number of unique genes in the dataset:", total_genes, "\n")

# Check the total number of rows in the global dataset (includes duplicates across samples/conditions)
total_rows <- nrow(res_global)
cat("Total rows in res_global (including duplicates):", total_rows, "\n")

# Double-check the number of unique genes for validation
unique_genes <- nrow(res_global %>% distinct(Gene))
cat("Number of unique genes confirmed in res_global:", unique_genes, "\n")











