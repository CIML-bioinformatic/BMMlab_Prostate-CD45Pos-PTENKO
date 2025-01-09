# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

## @knitr loadData

# Create the seurat object (RNA only)
merge_seurat_object = readRDS( PATH_SEURAT_OBJECT_INTEGRATION)

cat("<BR>Dispersion of cells among timepoints and conditions")
table( merge_seurat_object$orig.condition, merge_seurat_object$orig.date) %>% kable() %>% kable_styling( full_width = FALSE)
