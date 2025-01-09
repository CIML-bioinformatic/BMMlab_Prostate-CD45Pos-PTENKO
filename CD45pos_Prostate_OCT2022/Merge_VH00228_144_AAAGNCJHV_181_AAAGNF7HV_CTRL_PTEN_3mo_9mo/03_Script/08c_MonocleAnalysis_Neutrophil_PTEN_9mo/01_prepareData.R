# ################################################
# This script aims to read and filter data
# ################################################


## @knitr prepare_data

# READ THE DATA
# ...............

## @knitr prepare_data

# Load the Seurat object from heterogeneity analysis
cat("<BR><b>Loading data from\n", SEURAT_OBJECT_HETEROGENEITY_ANALYSIS, "</b><BR>")

sc10x.seurat = readRDS( SEURAT_OBJECT_HETEROGENEITY_ANALYSIS)

# Look at the cell numbers
print(
  table( sc10x.seurat$orig.condition, sc10x.seurat$orig.date) %>% kable() %>% kable_styling( full_width = FALSE)
)
cat("<BR>Total number of cells =", length(  Cells( sc10x.seurat)))

# Define the cells to select
selected_cells_condition = Cells( sc10x.seurat)[ sc10x.seurat[[ "orig.condition"]] == SELECTED_CONDITION]
selected_cells_date = Cells( sc10x.seurat)[ sc10x.seurat[[ "orig.date"]] == SELECTED_DATE]
selected_cells = intersect( selected_cells_condition, selected_cells_date)

cat("<BR><b>WARNING : Selecting cells from condition", SELECTED_CONDITION, "at date", SELECTED_DATE,"</b>")
cat("<BR>Number of cells of condition", SELECTED_CONDITION,"=", length( selected_cells_condition))
cat("<BR>Number of cells of date", SELECTED_DATE,"=", length( selected_cells_date))
cat("<BR>Number of selected cells =", length( selected_cells))
cat("<BR>Number of excluded cells =", length( Cells( sc10x.seurat)) - length( selected_cells))

# Filter the Seurat object to keep only the selected cells
sc10x.seurat = subset( sc10x.seurat, cells = selected_cells)

# Convert Seurat Object to cell_data_set object
cds_sc10X_seurat = as.cell_data_set( sc10x.seurat)

# Select the  module of genes the score of which was computed in the previous analysis steps
modified_module_names = paste0( names( MODULES_GENES), "1")
kept_modules = names( MODULES_GENES)[ which( modified_module_names %in% names( sc10x.seurat@meta.data))]
MODULES_GENES = MODULES_GENES[  kept_modules]
