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

if( QC_EXPLORATION_MODE == TRUE){
  outlier_vs_predicted.celltype.l1_table = table( merge_seurat_object_for_azimuth$outlier, merge_seurat_object_for_azimuth$predicted.celltype.l1)
  outlier_vs_predicted.celltype.l1_table %>% kable( caption = "Cell Outlier status versus cell type prediction L1") %>% kable_styling()
  
  outlier_vs_predicted.celltype.l1_chisq_test = chisq.test( outlier_vs_predicted.celltype.l1_table)
  corrplot::corrplot( outlier_vs_predicted.celltype.l1_chisq_test$residuals, is.cor = FALSE)
}



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

# ............................................................................
## Try to analyse the clusters to affect to each of them a dominant cell type
# ............................................................................

cat(" \n \n")
cat("### Prediction of cell type by cluster {.tabset .tabset-fade}")
cat(" \n \n")

# Compute and display the cluster prediction at cell type level 1
# .................................................................

cat(" \n \n")
cat("#### Prediction level 1 {.tabset .tabset-fade}")
cat(" \n \n")


# Assign to each cluster the cell type with the highest Pearson Residual
Idents( merge_seurat_object_for_azimuth) = "seurat_clusters"
cluster_identity = vector()
for( cluster in levels( Idents( merge_seurat_object_for_azimuth))){
  cluster_identity[ cluster] = as.character( names( which.max( predicted.celltype.l1_vs_cluster_chisq_test$residuals[ , cluster])))
}

kable( cluster_identity, caption = "Assigned L1 cell type by cluster", col.names = c("Cluster", "L1 cell type")) %>% kable_styling()
  
# Assign to each cell the cell type that has been assigned to its cluster
cell_cluster_identity = unlist( sapply( Idents( merge_seurat_object_for_azimuth), function( cluster){
  return( cluster_identity[ as.character( cluster)])
}))
names( cell_cluster_identity) = names(  Idents( merge_seurat_object_for_azimuth))
merge_seurat_object_for_azimuth = AddMetaData( merge_seurat_object_for_azimuth, metadata = cell_cluster_identity, col.name = "cell.type.cluster.identity.L1")

# Transfer the L1 cluster prediction to the original Seurat object
merge_seurat_object = AddMetaData( merge_seurat_object, 
                                   metadata = merge_seurat_object_for_azimuth$cell.type.cluster.identity.L1[ Cells( merge_seurat_object)],
                                   col.name = "cell.type.cluster.identity.L1") 

# Show the UMAP with cluster cell type identify
DimPlot( merge_seurat_object, group.by = "cell.type.cluster.identity.L1", label = TRUE, label.size = 3) + ggtitle( "Cells by cluster")

Idents( merge_seurat_object_for_azimuth) = "cell.type.cluster.identity.L1"
invisible( lapply( levels( Idents( merge_seurat_object_for_azimuth)), function( cell_type)
{
  
  cat(" \n \n")
  cat("##### ", cell_type)
  cat(" \n \n")
  
  print(
    tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
      DimPlot( merge_seurat_object, 
               reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
               cells.highlight = WhichCells( merge_seurat_object_for_azimuth, idents = cell_type),
               order = TRUE, 
               label = FALSE, 
               label.size = 6)  +
        ggtitle( paste0( cell_type, " (by dominant cell type in cluster)" )) + 
        theme( axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "none")),
      error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  
  cat(" \n \n"); # Required for '.tabset'
}))
cat(" \n \n"); # Required for '.tabset'

# Compute and display the cluster prediction at cell type level 2
# .................................................................

cat(" \n \n")
cat("#### Prediction level 2 {.tabset .tabset-fade}")
cat(" \n \n")

# Assign to each cluster the cell type with the highest Pearson Residual
Idents( merge_seurat_object_for_azimuth) = "seurat_clusters"
cluster_identity = vector()
for( cluster in levels( Idents( merge_seurat_object_for_azimuth))){
  cluster_identity[ cluster] = as.character( names( which.max( predicted.celltype.l2_vs_cluster_chisq_test$residuals[ , cluster])))
}

kable( cluster_identity, caption = "Assigned L2 cell type by cluster", col.names = c("Cluster", "L2 cell type")) %>% kable_styling()

# Assign to each cell the cell type that has been assigned to its cluster
cell_cluster_identity = unlist( sapply( Idents( merge_seurat_object_for_azimuth), function( cluster){
  return( cluster_identity[ as.character( cluster)])
}))
names( cell_cluster_identity) = names(  Idents( merge_seurat_object_for_azimuth))
merge_seurat_object_for_azimuth = AddMetaData( merge_seurat_object_for_azimuth, metadata = cell_cluster_identity, col.name = "cell.type.cluster.identity.L2")

# Transfer the L2 cluster prediction to the original Seurat object
merge_seurat_object = AddMetaData( merge_seurat_object, 
                                   metadata = merge_seurat_object_for_azimuth$cell.type.cluster.identity.L2[ Cells( merge_seurat_object)],
                                   col.name = "cell.type.cluster.identity.L2") 

# Show the UMAP with cluster cell type identify
DimPlot( merge_seurat_object, group.by = "cell.type.cluster.identity.L2", label = TRUE, label.size = 3) + ggtitle( "Cells by cluster")

Idents( merge_seurat_object_for_azimuth) = "cell.type.cluster.identity.L2"
invisible( lapply( levels( Idents( merge_seurat_object_for_azimuth)), function( cell_type)
{
  
  cat(" \n \n")
  cat("##### ", cell_type)
  cat(" \n \n")
  
  print(
    tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
      DimPlot( merge_seurat_object, 
               reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
               cells.highlight = WhichCells( merge_seurat_object_for_azimuth, idents = cell_type),
               order = TRUE, 
               label = FALSE, 
               label.size = 6)  +
        ggtitle( paste0( cell_type, " (by dominant cell type in cluster)" )) + 
        theme( axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "none")),
      error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  
  cat(" \n \n"); # Required for '.tabset'
}))
cat(" \n \n"); # Required for '.tabset'

# Compute and display the cluster prediction at cell type level 3
# .................................................................

cat(" \n \n")
cat("#### Prediction level 3 {.tabset .tabset-fade}")
cat(" \n \n")

# Assign to each cluster the cell type with the highest Pearson Residual
Idents( merge_seurat_object_for_azimuth) = "seurat_clusters"
cluster_identity = vector()
for( cluster in levels( Idents( merge_seurat_object_for_azimuth))){
  cluster_identity[ cluster] = as.character( names( which.max( predicted.celltype.l3_vs_cluster_chisq_test$residuals[ , cluster])))
}

kable( cluster_identity, caption = "Assigned L3 cell type by cluster", col.names = c("Cluster", "L3 cell type")) %>% kable_styling()

# Assign to each cell the cell type that has been assigned to its cluster
cell_cluster_identity = unlist( sapply( Idents( merge_seurat_object_for_azimuth), function( cluster){
  return( cluster_identity[ as.character( cluster)])
}))
names( cell_cluster_identity) = names(  Idents( merge_seurat_object_for_azimuth))
merge_seurat_object_for_azimuth = AddMetaData( merge_seurat_object_for_azimuth, metadata = cell_cluster_identity, col.name = "cell.type.cluster.identity.L3")

# Transfer the L3 cluster prediction to the original Seurat object
merge_seurat_object = AddMetaData( merge_seurat_object, 
                                   metadata = merge_seurat_object_for_azimuth$cell.type.cluster.identity.L3[ Cells( merge_seurat_object)],
                                   col.name = "cell.type.cluster.identity.L3") 

# Show the UMAP with cluster cell type identify
DimPlot( merge_seurat_object, group.by = "cell.type.cluster.identity.L3", label = TRUE, label.size = 3) + ggtitle( "Cells by cluster")

Idents( merge_seurat_object_for_azimuth) = "cell.type.cluster.identity.L3"
invisible( lapply( levels( Idents( merge_seurat_object_for_azimuth)), function( cell_type)
{
  
  cat(" \n \n")
  cat("##### ", cell_type)
  cat(" \n \n")
  
  print(
    tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
      DimPlot( merge_seurat_object, 
               reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
               cells.highlight = WhichCells( merge_seurat_object_for_azimuth, idents = cell_type),
               order = TRUE, 
               label = FALSE, 
               label.size = 6)  +
        ggtitle( paste0( cell_type, " (by dominant cell type in cluster)" )) + 
        theme( axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "none")),
      error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  
  cat(" \n \n"); # Required for '.tabset'
}))
cat(" \n \n"); # Required for '.tabset'


rm( "merge_seurat_object_for_azimuth")