#!/bin/bash
#SBATCH --job-name=run_fastp            # Job name
#SBATCH --partition=pibu_el8            # partition
#SBATCH --output=fastp_parallel_%j.out  # output log file
#SBATCH --error=fastp_parallel_%j.err   # error log file
#SBATCH --ntasks=1                      # Number of tasks 
#SBATCH --cpus-per-task=4               # Number of CPU cores to use
#SBATCH --time=04:00:00                 # Time limit for the job
#SBATCH --mem=32G                       # Total memory per job

# path to directory of the provided files 
READS_DIR="/data/courses/rnaseq_course/toxoplasma_de/reads"
# path to directory to store the resuls 
RESULTS_DIR="/data/users/sschaerer/rnaseq_course_files/fastp_results"
# container-path as found in the apptainer  
CONTAINER="/containers/apptainer/fastp_0.23.2--h5f740d0_3.sif"

# if non-existent, create an output directory for processed reads
mkdir -p "$RESULTS_DIR"

# Define a function to run Fastp on a single file
run_fastp() {
    # input defined as item following the command  
    file=$1
    # extracting the actual name from the input (removing th e extention)
    base=$(basename "$file" .fastq.gz)
    # defines the name of the output filepath by simply adding on _fastp to the extention
    output_file="${RESULTS_DIR}/${base}_fastp.fastq.gz"
    
    # --bind options map directories on the host system to directories inside the container
    apptainer exec --bind "$READS_DIR":"$READS_DIR" \
        --bind "$RESULTS_DIR":"$RESULTS_DIR" \
        "$CONTAINER" fastp \
        -i "$file" \                                    # input file 
        -o "$output_file" \                             # output file 
        -h "${RESULTS_DIR}/${base}_fastp.html" \        # HTML quality report
        -j "${RESULTS_DIR}/${base}_fastp.json"          # Json summary report 
}

# To minimise the runtime, the decision was made to run several instances of the function in parallel. 

# Export the function for GNU Parallel to be able to access it 
# GNU Parallel runs multiple instances of the function (run_fastp) in parallel
export -f run_fastp
export READS_DIR
export RESULTS_DIR
export CONTAINER

# find : generates a list of all .fastq.gz files in the READS_DIR
# | : pipes the list to GNU parallel for parallel processing 
# -j 4 : thefines the number of parallel running instances to 4 
# run_fastp : the previously defined function 
find "$READS_DIR"/*.fastq.gz | parallel -j 4 run_fastp
