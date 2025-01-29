
## STEP 5 : Exploratory data analysis

#-------------------------------------------------------------------------------

# script to create load the reformatted text file and create a DESeqDataSet object accordingly 
# runs the DESeq normalization
source("create_DESeqDataSet.R")


# Apply variance stabilizing transformation
# transforms the data to make the variance roughly constant across all expression levels
vst_data <- vst(dds, blind = TRUE)
# rows : genes 
# columns : samples 

# Generate PCA plot
# takes vst_data as input 
# uses group from sampleInfo as metadata to group and color the samples in the plot 
plotPCA(vst_data, intgroup = "group") +
  ggtitle("PCA of Grouped Samples") +
  theme(plot.title = element_text(hjust = 0.5))

# Notes for the plot 
# Points: Represent samples.
# Colors: Represent the group

# Each PC is a linear combination of the original variables (i.e. gene expression levels)
# -> each PC has an eigenvalue \lambda_k 
# percentage calculated = PC_k = \lambda_k/(sum of all \lambdas) * 100

# PC1/X-axis : captures the largest source of variation in the data.
# PC2 (y-axis) captures the second largest source of variation, orthogonal to PC1.


# Together : summary the major patterns of variability in the dataset.
# Percentages indicate how much of the total variability is explained by that component





