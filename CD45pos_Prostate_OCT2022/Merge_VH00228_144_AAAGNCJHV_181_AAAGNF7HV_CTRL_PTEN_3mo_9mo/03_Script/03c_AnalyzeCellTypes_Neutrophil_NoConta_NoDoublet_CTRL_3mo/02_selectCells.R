# #########################################
# This script aims to select only the cells
# of interest
# #########################################


# Select only the cells of interest
# --------------------------------------

## @knitr selectCells

Idents( merge_seurat_object) = "orig.condition"
merge_seurat_object = subset( merge_seurat_object, idents= CONDITION_NAME_CTRL)

Idents( merge_seurat_object) = "orig.date"
merge_seurat_object = subset( merge_seurat_object, idents= DATE_NAME_3MO)

cat("<BR>Dispersion of cells among timepoints and conditions after selection")
table( merge_seurat_object$orig.condition, merge_seurat_object$orig.date) %>% kable() %>% kable_styling( full_width = FALSE)

