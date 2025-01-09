# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

# CELL TYPE IDENTIFICATION BY AZIMUTH
#######################################

# ............................................
## @knitr heterogeneity_Azimuth
# ............................................

## Infer large cell types using Seurat Azimuth
# ............................................

cat("### Prediction of cell type on individual cells {.tabset .tabset-fade}")

# Convert mouse genes to human genes in count matrix
# ...................................................
count = GetAssayData( merge_seurat_object)
human_names_conversion = convertMouseGeneList( rownames( count))

human_name_set = vector()
for( gene_name in rownames( count)){
  index = which( human_names_conversion$symbol == gene_name)
  if( length( index) == 1){
    human_name_set = append( human_name_set, human_names_conversion[ index, "human_symbol"])
  }else if( length( index) > 1){
    human_name_set = append( human_name_set, human_names_conversion[ index[1], "human_symbol"])
  }else{
    human_name_set = append( human_name_set, NA)
  }
}

removed_rows = which( is.na( human_name_set))
count = count[ - removed_rows,]
human_name_set = human_name_set[ - removed_rows]

rownames( count) = human_name_set
merge_seurat_object_for_azimuth = CreateSeuratObject( counts = count)
merge_seurat_object_for_azimuth[[ "seurat_clusters"]] = merge_seurat_object[[ "seurat_clusters"]]

# Compute the Azimuth predictions
merge_seurat_object_for_azimuth <- RunAzimuth( merge_seurat_object_for_azimuth, reference = "pbmcref")

# Display the prediction at cell type level 1
# ............................................

cat(" \n \n")
cat("#### Prediction level 1")
cat(" \n \n")

# Transfer the L1 prediction to the original Seurat object
merge_seurat_object = AddMetaData( merge_seurat_object, 
                                   metadata = merge_seurat_object_for_azimuth$predicted.celltype.l1[ Cells( merge_seurat_object)],
                                   col.name = "predicted.celltype.l1") 

# Show the clusters on UMAP for reference
DimPlot( merge_seurat_object, group.by = "seurat_clusters", label = TRUE, label.size = 3) + ggtitle( "Cell clusters")

# Show the Azimuth result on UMAP
DimPlot( merge_seurat_object, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3) + ggtitle( "Azimuth PBMC Level 1 cell type assignation")
cat(" \n \n")

# Analyse the distribution of cell type against clusters
predicted.celltype.l1_vs_cluster_table = table( merge_seurat_object_for_azimuth$predicted.celltype.l1, merge_seurat_object_for_azimuth$seurat_clusters)
predicted.celltype.l1_vs_cluster_table %>% kable( caption = "Cell distribution of Azimuth prediction L1 versus Clusters") %>% kable_styling()

predicted.celltype.l1_vs_cluster_chisq_test = chisq.test( predicted.celltype.l1_vs_cluster_table)
corrplot::corrplot( predicted.celltype.l1_vs_cluster_chisq_test$residuals, is.cor = FALSE)


cat(" \n \n"); # Required for '.tabset'

# Display the prediction at cell type level 2
# ............................................

cat(" \n \n")
cat("#### Prediction level 2")
cat(" \n \n")

# Transfer the L2 prediction to the original Seurat object
merge_seurat_object = AddMetaData( merge_seurat_object, 
                                   metadata = merge_seurat_object_for_azimuth$predicted.celltype.l2[ Cells( merge_seurat_object)],
                                   col.name = "predicted.celltype.l2") 

# Show the clusters on UMAP for reference
DimPlot( merge_seurat_object, group.by = "seurat_clusters", label = TRUE, label.size = 3) + ggtitle( "Cell clusters")

# Show the Azimuth result on UMAP
DimPlot( merge_seurat_object, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + ggtitle( "Azimuth PBMC Level 2 cell type assignation")
cat(" \n \n")

# Analyse the distribution of cell type against clusters
predicted.celltype.l2_vs_cluster_table = table( merge_seurat_object_for_azimuth$predicted.celltype.l2, merge_seurat_object_for_azimuth$seurat_clusters)
predicted.celltype.l2_vs_cluster_table %>% kable( caption = "Cell distribution of Azimuth prediction L2 versus Clusters") %>% kable_styling()

predicted.celltype.l2_vs_cluster_chisq_test = chisq.test( predicted.celltype.l2_vs_cluster_table)
corrplot::corrplot( predicted.celltype.l2_vs_cluster_chisq_test$residuals, is.cor = FALSE)


cat(" \n \n"); # Required for '.tabset'

# Display the prediction at cell type level 3
# ............................................

cat(" \n \n")
cat("#### Prediction level 3")
cat(" \n \n")

# Transfer the L3 prediction to the original Seurat object
merge_seurat_object = AddMetaData( merge_seurat_object, 
                                   metadata = merge_seurat_object_for_azimuth$predicted.celltype.l3[ Cells( merge_seurat_object)],
                                   col.name = "predicted.celltype.l3") 

# Show the clusters on UMAP for reference
DimPlot( merge_seurat_object, group.by = "seurat_clusters", label = TRUE, label.size = 3) + ggtitle( "Cell clusters")

# Show the Azimuth result on UMAP
DimPlot( merge_seurat_object, group.by = "predicted.celltype.l3", label = TRUE, label.size = 3) + ggtitle( "Azimuth PBMC Level 3 cell type assignation")
cat(" \n \n")

# Analyse the distribution of cell type against clusters
predicted.celltype.l3_vs_cluster_table = table( merge_seurat_object_for_azimuth$predicted.celltype.l3, merge_seurat_object_for_azimuth$seurat_clusters)
predicted.celltype.l3_vs_cluster_table %>% kable( caption = "Cell distribution of Azimuth prediction L3 versus Clusters") %>% kable_styling()

predicted.celltype.l3_vs_cluster_chisq_test = chisq.test( predicted.celltype.l3_vs_cluster_table)
corrplot::corrplot( predicted.celltype.l3_vs_cluster_chisq_test$residuals, is.cor = FALSE)


cat(" \n \n"); # Required for '.tabset'


rm( "merge_seurat_object_for_azimuth")