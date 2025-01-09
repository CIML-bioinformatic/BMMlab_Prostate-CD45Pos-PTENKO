# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

## @knitr loadData

# Create the seurat object (RNA only)

cat("#### ANALYZING CELL TYPE :", SELECTED_CELL_TYPE)
cat("\n")
cat("\nLoad data from:\n", PATH_SEURAT_OBJECT_RDS_FILE)

# Load Seurat object from file
merge_seurat_object = readRDS( PATH_SEURAT_OBJECT_RDS_FILE)

  
