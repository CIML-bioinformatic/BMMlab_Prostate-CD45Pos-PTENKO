# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

# .....................................................
## @knitr loadData
# .....................................................

# Create the seurat object (RNA only)
merge_seurat_object = readRDS( file = PATH_SEURAT_OBJECT_RDS_FILE)

