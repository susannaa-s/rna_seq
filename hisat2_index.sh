#!/bin/bash
#SBATCH --job-name=hisat2_index         # Job name
#SBATCH --partition=pibu_el8            # Specify the partition
#SBATCH --output=hisat2_index_%j.out    # Standard output log file
#SBATCH --error=hisat2_index_%j.err     # Standard error log file
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPU cores to use
#SBATCH --time=04:00:00                 # Time limit for the job
#SBATCH --mem=32G                       # Total memory per job

# Define file paths and URLs
FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
# oath to output directory 
OUTPUT_DIR="/data/users/sschaerer/rnaseq_course_files/hisat2_index"
# defining a prefix to identify the file 
OUTPUT_PREFIX="$OUTPUT_DIR/Mus_musculus.GRCm39"
# defining the full path and file name for the decompressed version of the FASTA file
TEMP_FASTA="$OUTPUT_DIR/$(basename "$FASTA_URL" .gz)"
# path to container-tool as provided 
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Download the FASTA file if it doesn't already exist
# this has proven to be produce better results than downloading the file previously
# the resaons for this are unknown 
if [[ ! -f "$TEMP_FASTA.gz" ]]; then
    wget -O "$TEMP_FASTA.gz" "$FASTA_URL"
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
else
    echo "FASTA file already exists: $TEMP_FASTA.gz"
fi

# Decompress the FASTA file 
echo "Decompressing FASTA file to: $TEMP_FASTA"
gunzip -c "$TEMP_FASTA.gz" > "$TEMP_FASTA"

# Run hisat2-build using the decompressed FASTA file to create the genome index
echo "Running Hisat2 index creation..."
apptainer exec \
    --bind /data/users/sschaerer/rnaseq_course_files:/data/users/sschaerer/rnaseq_course_files \
    "$CONTAINER" hisat2-build -p 8 "$TEMP_FASTA" "$OUTPUT_PREFIX"

# Check for errors during index creation
if [[ $? -ne 0 ]]; then
    echo "Error: Hisat2 index creation failed."
    exit 1
fi

echo "Hisat2 index files have been successfully created in: $OUTPUT_DIR"
