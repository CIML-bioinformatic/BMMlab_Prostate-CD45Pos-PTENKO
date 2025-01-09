# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

## @knitr loadData

# Create the seurat object (RNA only)

sc10x_seurat_list = list()
for( sample_name in SAMPLE_NAME_LIST){
  
  sc10x_seurat_list[[ sample_name]] = CreateSeuratObject( counts = Read10X( PATH_SAMPLE_10X_FILES_LIST[[ sample_name]]),
                                                                  min.cells = LOAD_MIN_CELLS,
                                                                  min.features = LOAD_MIN_FEATURES,
                                                                  project = sample_name);

}

