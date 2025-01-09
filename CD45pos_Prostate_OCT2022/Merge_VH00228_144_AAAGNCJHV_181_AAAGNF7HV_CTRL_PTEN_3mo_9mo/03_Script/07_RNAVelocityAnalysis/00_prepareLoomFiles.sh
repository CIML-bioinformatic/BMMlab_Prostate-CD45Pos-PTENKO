#!/bin/bash

# This script aims to execute the velocyto program.
# 

# --------------------------------------------------
# Merge the BAM files from the various experiments
# --------------------------------------------------

# Define the location of the BAM files
BAM_CTRL3mo_1="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_1/05_Output/01_CellRangerAnalysis/mm10/outs/possorted_genome_bam.bam"
BAM_CTRL3mo_2="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_2/05_Output/01_CellRangerAnalysis/mm10/outs/possorted_genome_bam.bam"
BAM_PTEN3mo_1="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN3mo_1/05_Output/01_CellRangerAnalysis/mm10/outs/possorted_genome_bam.bam"
BAM_PTEN3mo_2="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN3mo_2/05_Output/01_CellRangerAnalysis/mm10/outs/possorted_genome_bam.bam"
BAM_CTRL9mo_1="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL9mo_1/05_Output/01_CellRangerAnalysis/mm10/outs/possorted_genome_bam.bam"
BAM_CTRL9mo_2="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL9mo_2/05_Output/01_CellRangerAnalysis/mm10/outs/possorted_genome_bam.bam"
BAM_PTEN9mo_1="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN9mo_1/05_Output/01_CellRangerAnalysis/mm10/outs/possorted_genome_bam.bam"
BAM_PTEN9mo_2="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN9mo_2/05_Output/01_CellRangerAnalysis/mm10/outs/possorted_genome_bam.bam"

# Define the location of the barcodes files
BARCODES_CTRL3mo_1="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_1/05_Output/01_CellRangerAnalysis/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BARCODES_CTRL3mo_2="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_2/05_Output/01_CellRangerAnalysis/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BARCODES_PTEN3mo_1="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN3mo_1/05_Output/01_CellRangerAnalysis/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BARCODES_PTEN3mo_2="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN3mo_2/05_Output/01_CellRangerAnalysis/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BARCODES_CTRL9mo_1="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL9mo_1/05_Output/01_CellRangerAnalysis/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BARCODES_CTRL9mo_2="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL9mo_2/05_Output/01_CellRangerAnalysis/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BARCODES_PTEN9mo_1="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN9mo_1/05_Output/01_CellRangerAnalysis/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BARCODES_PTEN9mo_2="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/VH00228_144_AAAGNCJHV_181_AAAGNF7HV_PTEN9mo_2/05_Output/01_CellRangerAnalysis/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

# Define the location of the output file
OUTPUT_FOLDER="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo/05_Output/07_RNAVelocityAnalysis/VelocytoLoom"

# Define the name of the merged BAM file
MERGED_BAM_FILE_NAME="Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo.bam"
MERGED_SORTED_BAM_FILE_NAME="Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo_sorted.bam"

# Define the name of the merged barcodes file
MERGED_BARCODES_FILE_NAME="Merge_Barcodes_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo.tsv"

# Provide the path to the genome annotation file
GENOME_ANNOTATION_PATH="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/01_REFERENCE/01_GENOME/CellRanger/mouse/2020-A/refdata-gex-mm10-2020-A/genes/genes.gtf"

# Create the output folder
mkdir -p $OUTPUT_FOLDER

cd $OUTPUT_FOLDER

# Merge the BAM in a single BAM files
mergeBams -i $BAM_CTRL3mo_1,$BAM_CTRL3mo_2,$BAM_PTEN3mo_1,$BAM_PTEN3mo_2,$BAM_CTRL9mo_1,$BAM_CTRL9mo_2,$BAM_PTEN9mo_1,$BAM_PTEN9mo_2 \
          -l CTRL3mo_1_,CTRL3mo_2_,PTEN3mo_1_,PTEN3mo_2_,CTRL9mo_1_,CTRL9mo_2_,PTEN9mo_1_,PTEN9mo_2_ \
          -b $BARCODES_CTRL3mo_1,$BARCODES_CTRL3mo_2,$BARCODES_PTEN3mo_1,$BARCODES_PTEN3mo_2,$BARCODES_CTRL9mo_1,$BARCODES_CTRL9mo_2,$BARCODES_PTEN9mo_1,$BARCODES_PTEN9mo_2 \
          -o $OUTPUT_FOLDER

# Rename the name of the merged bam file
mv $OUTPUT_FOLDER/out.bam ${OUTPUT_FOLDER}/${MERGED_BAM_FILE_NAME}

# unzip the barcodes files and add the correct prefix to the barcodes
gunzip -c $BARCODES_CTRL3mo_1 > $OUTPUT_FOLDER/barcodes_CTRL3mo_1.tsv
sed -i -e 's/^/CTRL3mo_1_/' $OUTPUT_FOLDER/barcodes_CTRL3mo_1.tsv

gunzip -c $BARCODES_CTRL3mo_2 > $OUTPUT_FOLDER/barcodes_CTRL3mo_2.tsv
sed -i -e 's/^/CTRL3mo_2_/' $OUTPUT_FOLDER/barcodes_CTRL3mo_2.tsv

gunzip -c $BARCODES_PTEN3mo_1 > $OUTPUT_FOLDER/barcodes_PTEN3mo_1.tsv
sed -i -e 's/^/PTEN3mo_1_/' $OUTPUT_FOLDER/barcodes_PTEN3mo_1.tsv

gunzip -c $BARCODES_PTEN3mo_2 > $OUTPUT_FOLDER/barcodes_PTEN3mo_2.tsv
sed -i -e 's/^/PTEN3mo_2_/' $OUTPUT_FOLDER/barcodes_PTEN3mo_2.tsv

gunzip -c $BARCODES_CTRL9mo_1 > $OUTPUT_FOLDER/barcodes_CTRL9mo_1.tsv
sed -i -e 's/^/CTRL9mo_1_/' $OUTPUT_FOLDER/barcodes_CTRL9mo_1.tsv

gunzip -c $BARCODES_CTRL9mo_2 > $OUTPUT_FOLDER/barcodes_CTRL9mo_2.tsv
sed -i -e 's/^/CTRL9mo_2_/' $OUTPUT_FOLDER/barcodes_CTRL9mo_2.tsv

gunzip -c $BARCODES_PTEN9mo_1 > $OUTPUT_FOLDER/barcodes_PTEN9mo_1.tsv
sed -i -e 's/^/PTEN9mo_1_/' $OUTPUT_FOLDER/barcodes_PTEN9mo_1.tsv

gunzip -c $BARCODES_PTEN9mo_2 > $OUTPUT_FOLDER/barcodes_PTEN9mo_2.tsv
sed -i -e 's/^/PTEN9mo_2_/' $OUTPUT_FOLDER/barcodes_PTEN9mo_2.tsv

# Concatenate all the prefixed barcodes files to a single file
cat $OUTPUT_FOLDER/barcodes_CTRL3mo_1.tsv $OUTPUT_FOLDER/barcodes_CTRL3mo_2.tsv $OUTPUT_FOLDER/barcodes_PTEN3mo_1.tsv $OUTPUT_FOLDER/barcodes_PTEN3mo_2.tsv $OUTPUT_FOLDER/barcodes_CTRL9mo_1.tsv $OUTPUT_FOLDER/barcodes_CTRL9mo_2.tsv $OUTPUT_FOLDER/barcodes_PTEN9mo_1.tsv $OUTPUT_FOLDER/barcodes_PTEN9mo_2.tsv > ${OUTPUT_FOLDER}/${MERGED_BARCODES_FILE_NAME}

# Ensure Samtools is in the PATH
export PATH=/samtools/bin:$PATH

# Sort the merged BAM file
samtools sort ${OUTPUT_FOLDER}/${MERGED_BAM_FILE_NAME} > ${OUTPUT_FOLDER}/${MERGED_SORTED_BAM_FILE_NAME}

# # Execute Velocyto on 10x data
echo "Command line is: velocyto run -b ${OUTPUT_FOLDER}/${MERGED_BARCODES_FILE_NAME} -o $OUTPUT_FOLDER ${OUTPUT_FOLDER}/${MERGED_BAM_FILE_NAME} $GENOME_ANNOTATION_PATH"
velocyto run -b ${OUTPUT_FOLDER}/${MERGED_BARCODES_FILE_NAME} -o $OUTPUT_FOLDER ${OUTPUT_FOLDER}/${MERGED_SORTED_BAM_FILE_NAME} $GENOME_ANNOTATION_PATH






