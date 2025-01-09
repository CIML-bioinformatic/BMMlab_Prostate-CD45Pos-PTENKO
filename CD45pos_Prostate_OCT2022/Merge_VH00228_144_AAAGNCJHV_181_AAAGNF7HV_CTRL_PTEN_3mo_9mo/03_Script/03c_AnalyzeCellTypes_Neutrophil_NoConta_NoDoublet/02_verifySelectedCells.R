# #########################################
# This script look at the selected cells
# to verify if they are the right one
# #########################################




# READ DATA
# ---------

## @knitr verifySelectedCells

# Check if the selected cells are the right one

Idents( merge_seurat_object) = "seurat_clusters"

table( Idents( merge_seurat_object)) %>% 
  kable( caption = "Selected cluster number and number of cells per clusters",col.names = c( "Cluster", "Nb of cells")) %>% 
  kable_styling( full_width = FALSE)

table(  merge_seurat_object$predicted.celltype.l1, Idents( merge_seurat_object)) %>% 
  kable( caption = "Selected cluster number and L1 predictited cell type (by cell)") %>% 
  kable_styling( full_width = FALSE)

table( merge_seurat_object$cell.type.cluster.identity.L1, Idents( merge_seurat_object)) %>% 
  kable( caption = "Selected cluster number and L1 predictited cell type (by cluster)") %>% 
  kable_styling( full_width = FALSE)
