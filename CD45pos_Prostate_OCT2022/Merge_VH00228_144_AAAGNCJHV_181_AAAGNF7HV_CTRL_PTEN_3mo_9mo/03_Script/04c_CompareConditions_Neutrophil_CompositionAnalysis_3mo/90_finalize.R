# ########################################################
# This script wraps up analysis (clean, save, render, ...)
# ########################################################




## @knitr final_saveSessionImage

# Save a binary file of final Seurat object only (as RDS)
seuratObjectPath = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "seuratObject_final.RDS"))
saveRDS( object = merge_seurat_object, file = seuratObjectPath)



