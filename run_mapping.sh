#!/bin/bash
#SBATCH --job-name=hisat2_mapping        # Job name
#SBATCH --partition=pibu_el8            # Partition
#SBATCH --output=hisat2_mapping_%j.out  # Standard output log file (%j inserts job ID)
#SBATCH --error=hisat2_mapping_%j.err   # Standard error log file
#SBATCH --ntasks=1                      # Number of tasks (single task)
#SBATCH --cpus-per-task=16              # Use 16 CPU cores for parallel processing
#SBATCH --time=08:00:00                 # Maximum runtime for the job
#SBATCH --mem=64G                       # Total memory allocated for the job

# Define paths for input, output, and reference files

# the decision was made to map the raw files, since : 
# 1. The quality check proved them to be good enough and 
# 2. The mapping with the trested files did not work due to differing lengths between the pairs 
# This step can be adjusted accordingly if the fastp treatment does not result in this problem
INPUT_DIR="/data/courses/rnaseq_course/toxoplasma_de/reads"                     # Directory containing raw paired-end reads

 # Directory to store mapped reads (BAM files)
OUTPUT_DIR="/data/users/sschaerer/rnaseq_course_files/mapped_reads"     
# Prefix for HISAT2 index files       
INDEX_PREFIX="/data/users/sschaerer/rnaseq_course_files/hisat2_index/Mus_musculus.GRCm39"  
# HISAT2 and Samtools container
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"          
# URL for GTF annotation file
GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz"  
# Path for the downloaded GTF file
GTF_FILE="/data/users/sschaerer/rnaseq_course_files/Mus_musculus.GRCm39.113.gtf.gz"      

# if non-existent, create an output directory
mkdir -p "$OUTPUT_DIR"

# Download the GTF file if it is not already present
# The GTF file is used for downstream analyses such as transcript quantification
if [[ ! -f "$GTF_FILE" ]]; then
    echo "Downloading GTF file from: $GTF_URL"
    wget -O "$GTF_FILE" "$GTF_URL"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to download GTF file."
        exit 1
    fi
else
    echo "GTF file already exists: $GTF_FILE"
fi

# Loop through all forward read files (_1.fastq.gz) in the input directory
for READ1 in "$INPUT_DIR"/*_1.fastq.gz; do
    READ2="${READ1/_1.fastq.gz/_2.fastq.gz}"  # Infer the reverse read file name (_2.fastq.gz)
    SAMPLE=$(basename "$READ1" _1.fastq.gz)   # Extract the sample name (base name without extensions)
    BAM_OUTPUT="$OUTPUT_DIR/${SAMPLE}.bam"   # Define the output BAM file path

    echo "Mapping reads for sample: $SAMPLE"

    # Perform the mapping process with HISAT2, followed by conversion to BAM and sorting using Samtools
    # HISAT2 options:
    # -p 16: Use 16 threads for parallel processing
    # --dta: Enable downstream transcript assembly
    # --rna-strandness RF: Specify stranded paired-end reads
    # Samtools:
    # - Convert SAM output from HISAT2 to BAM format
    # - Sort the BAM file by genomic coordinates for downstream compatibility
    apptainer exec \
        --bind /data/courses:/data/courses \
        --bind /data/users:/data/users \
        "$CONTAINER" hisat2 \
        -p 16 \
        --dta \
        --rna-strandness RF \
        -x "$INDEX_PREFIX" \
        -1 "$READ1" \
        -2 "$READ2" | \
    apptainer exec \
        --bind /data/courses:/data/courses \
        --bind /data/users:/data/users \
        "$CONTAINER" samtools view -@ 16 -bS - | \
    apptainer exec \
        --bind /data/courses:/data/courses \
        --bind /data/users:/data/users \
        "$CONTAINER" samtools sort -@ 16 -o "$BAM_OUTPUT"

    # Check for errors in the mapping process
    if [[ $? -ne 0 ]]; then
        echo "Error: Mapping failed for sample: $SAMPLE"
        exit 1
    fi

    echo "Mapping completed for sample: $SAMPLE"
done

