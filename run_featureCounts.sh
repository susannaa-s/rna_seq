#!/bin/bash
#SBATCH --job-name=featureCounts          # Job name
#SBATCH --partition=pibu_el8              # Partition to run the job
#SBATCH --output=featureCounts_%j.out     # Standard output log file (%j inserts the job ID)
#SBATCH --error=featureCounts_%j.err      # Standard error log file
#SBATCH --ntasks=1                        # Number of tasks (single task)
#SBATCH --cpus-per-task=8                 # Number of CPU cores allocated
#SBATCH --time=04:00:00                   # Maximum runtime for the job
#SBATCH --mem=32G                         # Total memory allocated for the job

# Define directories and file paths

# Directory containing sorted, indexed BAM files
BAM_DIR="/data/users/sschaerer/rnaseq_course_files/s_i_m_reads" 
# Path to the GTF annotation file
ANNOTATION_FILE="/data/users/sschaerer/rnaseq_course_files/Mus_musculus.GRCm39.113.gtf.gz"  
# Path for the output count matrix
OUTPUT_FILE="/data/users/sschaerer/rnaseq_course_files/featureCounts_counts.txt"  
# Path to the container with FeatureCounts
CONTAINER="/containers/apptainer/subread_2.0.1--hed695b0_0.sif"         

# Gather all BAM files in the directory into a single variable
# 'find' locates all BAM files, and 'tr' converts newline-separated paths into space-separated paths
BAM_FILES=$(find "$BAM_DIR" -name "*.bam" | tr '\n' ' ')

# Ensure proper read permissions for the BAM files
# chmod +r ensures read permissions for BAM files, and chmod +rx ensures directory permissions
chmod +r "$BAM_DIR"/*.bam
chmod +rx "$BAM_DIR"

# Run FeatureCounts to assign reads to genes
# Options:
# -T 8: Use 8 threads for parallel processing
# -a: Specify the GTF annotation file for gene assignment
# -o: Define the output file for the read count matrix
# -s 2: Specify reverse-strandedness of the data
# -p: Enable paired-end read support
apptainer exec \
    --bind /data/users/sschaerer/rnaseq_course_files:/data/users/sschaerer/rnaseq_course_files \
    "$CONTAINER" \
    featureCounts -T 8 -a "$ANNOTATION_FILE" -o "$OUTPUT_FILE" -s 2 -p $BAM_FILES

# Confirm completion of read counting
echo "Read counting completed. Output written to $OUTPUT_FILE"



