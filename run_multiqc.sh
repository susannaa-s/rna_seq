#!/bin/bash
#SBATCH --job-name=run_multiqc          # Job name
#SBATCH --partition=pibu_el8            # Partition where the job will run
#SBATCH --output=multiqc_%j.out         # Output log file (%j inserts job ID)
#SBATCH --error=multiqc_%j.err          # Error log file
#SBATCH --ntasks=1                      # Number of tasks (one task for this job)
#SBATCH --cpus-per-task=4               # Number of CPU cores allocated
#SBATCH --time=01:00:00                 # Maximum runtime for the job
#SBATCH --mem=16G                       # Total memory allocated for the job

# To make this script general and reusable, input and output directory paths are specified as $1 (input) and $2 (output).
INPUT_DIR=$1
OUTPUT_DIR=$2

# if non-existent, create an output directory
mkdir -p "$OUTPUT_DIR"

# Path to the MultiQC container
MULTIQC_CONTAINER="/containers/apptainer/multiqc-1.19.sif"

# Run MultiQC using the specified container
echo "Running MultiQC on files in $INPUT_DIR..."
apptainer exec --bind "$INPUT_DIR":"$INPUT_DIR","$OUTPUT_DIR":"$OUTPUT_DIR" \ # Map input and output directories into the container
    "$MULTIQC_CONTAINER" multiqc \                                           # Run MultiQC inside the container
    "$INPUT_DIR" -o "$OUTPUT_DIR"                                            # Specify the input directory and output directory for MultiQC

# Running MultiQC is a relatively quick process, so it was not necessary to define a function or run multiple instances in parallel with GNU Parallel.


