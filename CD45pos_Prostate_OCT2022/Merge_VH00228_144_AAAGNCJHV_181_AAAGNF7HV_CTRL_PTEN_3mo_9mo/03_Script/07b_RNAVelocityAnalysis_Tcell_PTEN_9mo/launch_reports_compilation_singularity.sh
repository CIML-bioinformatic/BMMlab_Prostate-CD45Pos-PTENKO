#!/bin/bash



# The path to the Singularity Image
SINGULARITY_IMAGE_PATH="/mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo/02_Container/bmmlab_prostatebmshdp_velocyto.sif"

# The path to the script to execute
SCRIPT_FOLDER_PATH="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/spinelli/ciml-bmmlab/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo/03_Script/07_RNAVelocityAnalysis/"

# Change dir to the script folder
cd $SCRIPT_FOLDER_PATH

# Launch the script through in a Singularity image
singularity exec -B /mnt:/mnt $SINGULARITY_IMAGE_PATH Rscript ./launch_reports_compilation.R