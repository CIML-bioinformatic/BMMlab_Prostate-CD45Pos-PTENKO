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

## Remove doublet cells
Idents( merge_seurat_object) = "scDblFinder"
doublet_cells = Cells( merge_seurat_object)[ which( Idents( merge_seurat_object) == "doublet")]
write.table( doublet_cells, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, "Doublet_cells.csv"))
  
# Remove the doublet cells
cat("<BR><b>Removing", length( doublet_cells), "doublet cells</b><BR>")
merge_seurat_object = subset( merge_seurat_object, cells = setdiff( Cells( merge_seurat_object), doublet_cells))
  
