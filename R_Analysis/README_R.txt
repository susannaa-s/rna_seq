Differential Gene Expression Analysis Workflow

This README provides a step-by-step guide for reproducing the differential expression analysis workflow, starting from raw featureCounts output. The featureCount files as well as all the resulting files and scripts should be stored in the same directory, serving as the working directory when editing them with R. 

For detailed explanations of each step, please refer to the comments in the corresponding scripts. For a general understanding of the process and analysis of the results, consult the corresponding section in the project report PDF.

---

Before anything please unzipp the following zipped files in the directory. (zipping was done due to size limit on github) 
- featureCounts_counts.txt.zip 
- remormatted_counts.txt.zip

Step 1: Reformatting the featureCounts_counts.txt File
Objective:
Prepare the featureCounts_counts.txt file to match the input structure expected by DESeq2.

Procedure:
1. Navigate to the directory containing the raw featureCounts_counts.txt file:
   
   cd "/Users/susannascharer/Desktop/Master Semester 1/RNA sequencing/Group Project"
   
2. Remove the header and unnecessary columns (`Chr`, `Start`, `End`, `Strand`, and `Length`) using the following command:
   
   tail -n +2 "featureCounts_counts.txt" | cut --complement -f 2-6 > "R_Analysis/reformatted_counts.txt"
   
3. The reformatted counts matrix is now saved as `reformatted_counts.txt` in the R_Analysis directory.

---

Step 2: Create a DESeqDataSet Object
Objective:
Generate a DESeqDataSet object using the reformatted counts and associated sample metadata. Then perform differential expression analysis to identify genes with significant expression changes across conditions.

Requirements:
- Reformatted counts file: `reformatted_counts.txt`
- Information on the sample metadata

This file will be run within every subsequent file that is created. Therefore, the decision was made to generalise these steps into a specific script.

- Script Name: create_DESeqDataSet.R
- Full Path: ~/Desktop/Master Semester 1/RNA sequencing/Group Project/R_Analysis/create_DESeqDataSet.R

At this point, there is no need to run the script directly.

---

Step 3: Quality Control and Visualisation
Objective:
Inspect and visualise the DESeq2 results.

Since the `create_DESeqDataSet.R` script is sourced in the first line of this step, it is recommended to inspect the resulting DESeqDataSet object. Visualise principal component analysis (PCA) by running the following script:

- Script Name: rnaseq_Step5.R
- Full Path: ~/Desktop/Master Semester 1/RNA sequencing/Group Project/R_Analysis/rnaseq_Step5.R
  

Step 4: Differential Expression Analysis
Objective: Identify significant genes across conditions and rank them based on adjusted p-values ($p_{adj}$) and log$_2$ fold changes. Analyse genes globally and within specific tissue datasets (blood and lung).

Script :
- Name: `rnaseq_step6.R`
- Path: `~/Desktop/Master Semester 1/RNA sequencing/Group Project/R_Analysis/rnaseq_step6.R`

Procedure:
1. Differential expression analysis was conducted using DESeq2 for the following comparisons:
   - Blood\_Case vs Blood\_Control
   - Lung\_Case vs Lung\_Control
   - Combined global analysis across all conditions
2. Genes were ranked and analysed based on:
   - Adjusted p-values ($p_{adj}$): For identifying statistically significant genes.
   - Absolute log$_2$ fold change: For prioritising biologically impactful genes.
3. Condition-specific metadata was integrated to ensure clear distinctions between Blood\_Case, Blood\_Control, Lung\_Case, and Lung\_Control.
4. Six genes of interest (\textit{Ifit1, Mx1, Oas1a, Gbp2, Tap1, Cd274}) identified from previous studies were evaluated across metrics for their statistical and biological significance.

Outputs (optional) :
- CSV files summarising:
  - Top 20 genes globally, for blood, and for lung based on $p_{adj}$ and log$_2$ fold change.
  - Rankings of the six selected genes globally and within blood and lung datasets.
- Volcano plots illustrating statistical and biological significance of genes across conditions.


Sources : 

Singhania, Akul et al. (2019). “Transcriptional profiling unveils type I and II interferon networks
in blood and tissues across diseases”. In: Nature Communications 10.1, p. 2887. doi: 10.1038/
s41467-019-10601-6. url: https://doi.org/10.1038/s41467-019-10601-6.








