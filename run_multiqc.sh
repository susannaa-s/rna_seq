#!/bin/bash
#SBATCH --job-name=run_multiqc
#SBATCH --partition=pibu_el8
#SBATCH --output=multiqc_%j.out
#SBATCH --error=multiqc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G

# directory with the fastqc files 
fastqc_results_dir="/data/users/${USER}/fastqc_results"
# directory to store the multiqc resulting files in 
multiqc_output_dir="/data/users/${USER}/multiqc_results"

# creating said dierectory if necessary 
mkdir -p "$multiqc_output_dir"

# we load the module if avaliable direcctly, otherwise we use the one avaliable in the appatainer container 
if module avail multiqc 2>/dev/null; then
    module load multiqc
    multiqc -o "$multiqc_output_dir" "$fastqc_results_dir"
else
    apptainer exec /data/users/${USER}/containers/multiqc-1.12.sif \
        multiqc -o "$multiqc_output_dir" "$fastqc_results_dir"
fi
