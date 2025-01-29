#!/bin/bash
#SBATCH --job-name=run_fastqc            # Job name
#SBATCH --partition=pibu_el8             # Partition for the job
#SBATCH --output=fastqc_parallel_%j.out  # Output log file (%j is the job ID)
#SBATCH --error=fastqc_parallel_%j.err   # Error log file
#SBATCH --ntasks=1                       # Number of tasks (single task)
#SBATCH --cpus-per-task=4                # Number of CPU cores allocated
#SBATCH --time=02:00:00                  # Time limit for the job
#SBATCH --mem=16G                        # Memory allocation

# To make this script general and reusable, input and output directory paths are specified as $1 (input) and $2 (output).

# Path to the FastQC container as provided
FASTQC_CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"

# Input and output directory paths (taken from command-line arguments)
INPUT_DIR="$1"
OUTPUT_DIR="$2"

# if non-existent, create an output directory
mkdir -p "$OUTPUT_DIR"

# Define a function to run FastQC on a single file
# This mirrors the structure of the run_fastp script for consistency
run_fastqc() {
    local file=$1
    echo "Running FastQC on $file"
    apptainer exec --bind "$INPUT_DIR":"$INPUT_DIR" \
                   --bind "$OUTPUT_DIR":"$OUTPUT_DIR" \
                   "$FASTQC_CONTAINER" fastqc \
                   -o "$OUTPUT_DIR" \
                   "$file"
}

# Export the function so that GNU Parallel can access and execute it
export -f run_fastqc

# To ensure variables are available to run_fastqc when executed in parallel, export them as well
export FASTQC_CONTAINER INPUT_DIR OUTPUT_DIR

# Find all .fastq.gz files in the input directory
# Pipe the file list to GNU Parallel for processing
# Limit the number of concurrent FastQC instances to 4 for optimal resource usage
find "$INPUT_DIR" -name "*.fastq.gz" | parallel -j 4 run_fastqc
