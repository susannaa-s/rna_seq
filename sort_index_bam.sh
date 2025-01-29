#!/bin/bash
#SBATCH --job-name=parallel_sort_bam    # Job name
#SBATCH --partition=pibu_el8            # Partition for the job
#SBATCH --output=parallel_sort_bam_%j.out  # Standard output log file (%j inserts job ID)
#SBATCH --error=parallel_sort_bam_%j.err   # Standard error log file
#SBATCH --ntasks=1                      # Number of tasks (single task)
#SBATCH --cpus-per-task=8               # Number of CPU cores allocated
#SBATCH --time=08:00:00                 # Maximum runtime for the job
#SBATCH --mem=32G                       # Total memory allocated for the job

# Define directories for input BAM files and output sorted BAM files

# Input directory containing unsorted BAM files
BAM_DIR="/data/users/sschaerer/rnaseq_course_files/mapped_reads"            
# Output directory for sorted and indexed BAM files
SORTED_BAM_DIR="/data/users/sschaerer/rnaseq_course_files/s_i_m_reads"      
# Path to the container with Samtools
CONTAINER="/containers/apptainer/hisat2_samtools_v4.0.0-beta.sif"           


# if non-existent, create the directory 
mkdir -p "$SORTED_BAM_DIR"

# Define a function to sort and index BAM files
sort_and_index_bam() {
    local bam_file=$1  # Input BAM file
    local sample_name=$(basename "$bam_file" .bam)  # Extract the sample name (without extension)
    local sorted_bam="$SORTED_BAM_DIR/${sample_name}_sorted.bam"  # Define path for the sorted BAM file

    echo "Processing $sample_name"

    # Sort the BAM file by genomic coordinates
    # -@ 2: Use 2 threads for parallel sorting within Samtools
    apptainer exec \
        --bind /data/users/sschaerer/rnaseq_course_files:/data/users/sschaerer/rnaseq_course_files \
        "$CONTAINER" samtools sort -@ 2 -o "$sorted_bam" "$bam_file"

    # Index the sorted BAM file
    apptainer exec \
        --bind /data/users/sschaerer/rnaseq_course_files:/data/users/sschaerer/rnaseq_course_files \
        "$CONTAINER" samtools index "$sorted_bam"

    echo "Processed $sample_name: Sorted BAM and Index created."
}

# Export the function and required variables for GNU Parallel
# This ensures that GNU Parallel can access the sort_and_index_bam function and related variables
export -f sort_and_index_bam
export SORTED_BAM_DIR CONTAINER

# Find all unsorted BAM files in the input directory
# Use GNU Parallel to process files in parallel with up to 4 concurrent jobs
find "$BAM_DIR" -name "*.bam" | parallel -j 4 sort_and_index_bam

# Final message indicating completion of the script
echo "All BAM files sorted and indexed."
