#!/bin/bash

# This script execute the CellRanger analysis of the Prostate_Project_BMSH_DP project
# for one of the project libraries. Sample nale is the name of the analyzed library


# Name of the library to analyse
SAMPLE_NAME="CTRL9mo_1"

# Path to the CellRanger singularity image
SINGULARITY_IMAGE="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/02_CONTAINER/02_SINGULARITY/CellRanger/7_0_1/cellranger_701.sif"

# Folder where the CellRanger reference file of the genome are stored
FOLDER_TO_REFERENCE="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/01_REFERENCE/01_GENOME/CellRanger/mouse/2020-A/refdata-gex-mm10-2020-A"

# Folder containing the links to the fastq files
FASTQ_FOLDER="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/230127_VH00228_144_AAAGNCJHV_${SAMPLE_NAME}/00_RawData/01_mRNA"

# Output folder for the CellRanger output files
OUTPUT_FOLDER="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/230127_VH00228_144_AAAGNCJHV_${SAMPLE_NAME}/05_Output/01_CellRangerAnalysis"

# Execution of the analysis
echo SAMPLE_NAME = $SAMPLE_NAME
echo FOLDER_TO_REFERENCE = $FOLDER_TO_REFERENCE
echo FASTQ_FOLDER = $FASTQ_FOLDER
echo OUTPUT_FOLDER = $OUTPUT_FOLDER
                                   
mkdir -p $OUTPUT_FOLDER

cd $OUTPUT_FOLDER
echo -e "\n\nCommand:\nsingularity exec -B /mnt/DOSI:/mnt/DOSI $SINGULARITY_IMAGE cellranger count --id=mm10 --transcriptome=$FOLDER_TO_REFERENCE --fastqs=$FASTQ_FOLDER --sample=$SAMPLE_NAME\n\n"
singularity exec -B /mnt/DOSI:/mnt/DOSI $SINGULARITY_IMAGE cellranger count --id=mm10 --transcriptome=$FOLDER_TO_REFERENCE --fastqs=$FASTQ_FOLDER --sample=$SAMPLE_NAME 





