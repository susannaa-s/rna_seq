#!/bin/bash
#SBATCH --job-name=run_fastqc           # Job name
#SBATCH --partition=pibu_el8            # Specify the partition
#SBATCH --output=fastqc_parallel_%j.out # Standard output log file
#SBATCH --error=fastqc_parallel_%j.err  # Standard error log file
#SBATCH --ntasks=1                      # Number of tasks 
#SBATCH --cpus-per-task=4               # Number of CPU cores to use
#SBATCH --time=02:00:00                 # Time limit for the job
#SBATCH --mem=16G                       # Total memory per job

# Create an output directory for FastQC results if it doesn't exist
mkdir -p /data/users/sschaerer/rnaseq_course/fastqc_results

# Define a function to run FastQC on a single file
# $1 refers to the first argument passed to the function when called
run_fastqc() {
    file=$1
    # Print a message to track the file being processed
    echo "Running FastQC on $file" 
    apptainer exec --bind /data/courses/rnaseq_course/toxoplasma_de/reads:/data/courses/rnaseq_course/toxoplasma_de/reads \
        /containers/apptainer/fastqc-0.12.1.sif fastqc -o /data/users/sschaerer/rnaseq_course/fastqc_results "$file"
}

# --bind : explicitly maps the reads directory to ensure the container can access it

# Export the function to make it accessible to GNU Parallel
export -f run_fastqc

# Find all fastq.gz files and run FastQC in parallel
find /data/courses/rnaseq_course/toxoplasma_de/reads/*.fastq.gz | parallel -j 4 run_fastqc

