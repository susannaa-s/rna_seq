# RNA sequencing project Workflow 

This file gives an overview of the sequence in files and how they should be applied in order to obtain the results. the context and discussion to the results can be found in the report pdf and will not be discussed here. This serves as manual, to replicate the result-generating process later on. 

Container tools used : 

- FastQC : `/containers/apptainer/fastqc-0.12.1.sif`
- MultiQC : `/containers/apptainer/multiqc-1.19.sif`
- Hisat2 : `/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif`
- Samtools : `/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif`
- featureCounts : `/containers/apptainer/subread_2.0.1--hed695b0_0.sif`

All the scripts are stored in the following directory : 
Name : rnaseq_course 
Complete path : /data/users/sschaerer/rnaseq_course


# Step 1 : Setup 

Once access has been gained to the cluster, we install conect a previously created Git-Repository to a folder in our Workspace called 'rnaseq_course'. This folder holds all necessary scripts to process the given files. 
Name : rnaseq_course 
Complete path : /data/users/sschaerer/rnaseq_course

We create another directory, which remains disconnected from the Git-Repository to hold all the large result files, since pushing those onto git has proven to be nearly impossible. 
Name : rnaseq_course_files 
Complete path : /data/users/sschaerer/rnaseq_course_files

# Step 2 : Quality Control 

We are given the path to a directory 
Name : READS
Complete path : /data/courses/rnaseq_course/toxoplasma_de/reads
with a total of 32 files, 2 per sample. More detailed information is provided in the local README file 
Name : README 
Complete path : /data/courses/rnaseq_course/toxoplasma_de/README

We start off the running FastP on all of the files. 

Once the terminal is opened in the propper directory (rnaseq_course) the following commands are to be run: 

1. `sbatch run_fastp.sh`

The files should be generated in a separate file 
Name : fastp_results 
Complete path : /data/users/sschaerer/rnaseq_course/fastp_results

We move on to running FastQC on the original and treated files by running the following commands. We have modified the file to be able to define the path to the input- and output directory, since we will be running it on the original and the processed files. 
We do so by running the following commands : 

2. `sbatch run_fastqc.sh /data/courses/rnaseq_course/toxoplasma_de/reads /data/users/sschaerer/rnaseq_course_files/fastqc_original_results`

3. `sbatch run_fastqc.sh /data/users/sschaerer/rnaseq_course_files/fastp_results /data/users/sschaerer/rnaseq_course_files/fastqc_fastp_results`


One both set of fies ave been generated in : 
Name : fastqc_original_results 
Complete path : /data/users/sschaerer/rnaseq_course_files/fastqc_original_results
Name : fastqc_fastp_results 
Complete path : /data/users/sschaerer/rnaseq_course_files/fastqc_fastp_results

We run multiqc on both of the sets  individually with the following commands : 

4. `sbatch run_multiqc.sh /data/users/sschaerer/rnaseq_course_files/fastqc_original_results /data/users/sschaerer/rnaseq_course_files/multiqc_original_results`

5. `sbatch run_multiqc.sh /data/users/sschaerer/rnaseq_course_files/fastqc_fastp_results /data/users/sschaerer/rnaseq_course_files/multiqc_fastp_results`

# Step 3 : Map reads to the reference genome 

We import the two indicated folders to the 'rna_course_files' directory. thex can be found under the following links : 
- https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
- https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz

We store both of them under the following terms : 
Name : Mus_musculus.GRCm39.113.gtf
Complete path : /data/users/sschaerer/rnaseq_course_files/Mus_musculus.GRCm39.113.gtf
Name : Mus_musculus.GRCm39.dna.primary_assembly.fa
Complete path : /data/users/sschaerer/rnaseq_course_files/Mus_musculus.GRCm39.dna.primary_assembly.fa

To produce the genome index required by HISAT2, the FASTA file is used as input for the hisat2-build command. The index organises the genome into a searchable structure, enabling efficient and accurate alignment of RNA-Seq reads. The process uses 8 CPU threads for parallelisation to reduce runtime. The command to initiate the process is:

6. `sbatch hisat2_index.sh`
- generates 8 resulting files (one per chromosome), Mus_musculus.GRCm39.X.ht2 (X in [1,8])
- can be found in the following directory 
Name : hisat2_index 
Complete path : /data/users/sschaerer/rnaseq_course_files/hisat2_index


We move on to the mapping process. Where we map each of the sample file-sets onto the previously generated index file with read strandedness RF (reverse-forward). Each of those mappings generate a .sam file which is converted into .bam format to save space. 

7. `run_mapping.sh`
- generates sixteen .bam files one for each sample 
- can be found int he following dorectory 
Name : mapped_reads 
Complete path : /data/users/sschaerer/rnaseq_course_files/mapped_reads

From this point on, we use Samtools to sort the bam files by genomic coordinates to Index the coordinate sorted bam files again using bam files.
To do so we run the following command : 
8. `sbatch sort_index_bam.sh`
- generates sixteen sorted - indexed - mapped - reads 
- can be found in the following directory 
Name : s_i_m_reads 
Complete path : /data/users/sschaerer/rnaseq_course_files/s_i_m_reads


# Step 4. Count the number of reads per gene

To gain further insights into the readings, we use featurecounts to generate a detailed report on the assignments status of differents reads in each of the samples. 
We do this by running the following script which applied featurecounts on each of the files stored in the s_i_m_reads directory and stires them in a new one. 

9. `sbatch run_featureCounts.sh`
- generates two files containing first the detailed results and second a summary of said results. Both are saved in the files-directory. Temporary files are also stored there during the running of the script. 

Name : featureCounts_counts.txt
Complete path : /data/users/sschaerer/rnaseq_course_files/featureCounts_counts.txt

Name : featureCounts_counts.txt.summary
Complete path : /data/users/sschaerer/rnaseq_course_files/featureCounts_counts.txt.summary

We run the multiQC script once more on these files to have better visual representation of the results and store them in a directory. 

10. `sbatch run_multiqc.sh /data/users/sschaerer/rnaseq_course_files /data/users/sschaerer/rnaseq_course_files/fc_multiQC_results`

Name : fc_multiQC_results 
Complete path : /data/users/sschaerer/rnaseq_course_files/fc_multiQC_results

A detailed look at those results can be found in the report and will not be discussed here. 


# END OF THE CLUSTER PART # 

At this point we have generated all the necessary files to collect for rither analysis. This will be done on the local machine using R. 

