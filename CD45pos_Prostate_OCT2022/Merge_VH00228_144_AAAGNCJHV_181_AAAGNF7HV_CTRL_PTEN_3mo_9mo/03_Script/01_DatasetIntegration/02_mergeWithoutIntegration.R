# #########################################
# This script execute a simple dataset
# merge with no integration
# #########################################


# MERGE WITH NO INTEGRATION
# ----------------------------

## @knitr mergeWithoutIntegration

# Merge the dataset into a single Seurat Object
merge_seurat_object <- merge( x = sc10x_seurat_list[[ SAMPLE_NAME_LIST[ 1] ]],
                                      y = c( sc10x_seurat_list[[ SAMPLE_NAME_LIST[ 2]]],
                                             sc10x_seurat_list[[ SAMPLE_NAME_LIST[ 3]]],
                                             sc10x_seurat_list[[ SAMPLE_NAME_LIST[ 4]]],
                                             sc10x_seurat_list[[ SAMPLE_NAME_LIST[ 5]]],
                                             sc10x_seurat_list[[ SAMPLE_NAME_LIST[ 6]]],
                                             sc10x_seurat_list[[ SAMPLE_NAME_LIST[ 7]]],
                                             sc10x_seurat_list[[ SAMPLE_NAME_LIST[ 8]]]), 
                                      add.cell.ids = SAMPLE_NAME_LIST, 
                                      project = INTEGRATION_SAMPLE_NAME
                                      )

# Remove individual samples Seurat objects to spare memory
for( sample_name in SAMPLE_NAME_LIST){
  sc10x_seurat_list[[ sample_name]] = NULL
}
rm ("sc10x_seurat_list")
gc()

# Attribute a numeric ID to each cell (easier than barcode)
merge_seurat_object[["numID"]] = 1:length( Cells( merge_seurat_object));


# Create an annotation with the condition of origin of each cell
# ....................................................................

Idents( merge_seurat_object) = "orig.ident"
condition_of_origin = sapply( Idents( merge_seurat_object), function( sample_name){
   return( as.character( SAMPLE_CONDITION[ sample_name]))
})

merge_seurat_object = AddMetaData( merge_seurat_object, metadata = condition_of_origin, col.name = "orig.condition")

# Create an annotation with the date of origin of each cell
# ....................................................................

Idents( merge_seurat_object) = "orig.ident"
date_of_origin = sapply( Idents( merge_seurat_object), function( sample_name){
  return( as.character( SAMPLE_DATE[ sample_name]))
})

merge_seurat_object = AddMetaData( merge_seurat_object, metadata = date_of_origin, col.name = "orig.date")

# Create an annotation cobining condition and date of origin of each cell
# .......................................................................
merge_seurat_object = AddMetaData( merge_seurat_object, metadata = paste( condition_of_origin, date_of_origin, sep="_"), col.name = "orig.condition_date")

