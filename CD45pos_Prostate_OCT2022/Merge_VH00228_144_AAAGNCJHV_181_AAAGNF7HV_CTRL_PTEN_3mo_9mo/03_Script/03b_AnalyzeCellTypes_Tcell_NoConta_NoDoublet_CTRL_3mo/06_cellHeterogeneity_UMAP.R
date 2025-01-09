# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# DIMENSIONAL REDUCTION (TSNE/UMAP)
###################################
# ..........................................................................................................
## @knitr dimreduc_ggplot_covariables
# ..........................................................................................................

# Plot the cells sample of origin in UMAP
DimPlot( merge_seurat_object, reduction = "umap", group.by = "orig.ident") +
  ggtitle( "Map of cells by sample of origin")

# Plot the cells sample by clusters
clustersColor = CLUSTER_COLOR_PANEL[ 1:length( unique( merge_seurat_object$seurat_clusters))]
names( clustersColor) = sort(  unique( merge_seurat_object$seurat_clusters))

DimPlot( merge_seurat_object, reduction = "umap", group.by = "seurat_clusters",
         cols = clustersColor) + 
  ggtitle( "Map of cells by cluster")

# if QC EXPLORATION MODE is active, plot if the cells are outliers or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( merge_seurat_object, reduction = "umap",
                  group.by = "outlier") + 
                  ggtitle( "Map of cells by \nQC filtering status")
  )
}

# Plot the cells percentage of ribosomal genes in UMAP
FeaturePlot( merge_seurat_object, features = "percent.ribo") +
              scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
              ggtitle( "Map of cells with level of\n percentage of ribosomal genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of ribosomal genes or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                                  group.by = "outlier.percent.ribo") + 
                 ggtitle( "Map of cells by filter status\n on %ribosomal gene")
      )
}

# Plot the cells percentage of mitochondrial genes in UMAP  
FeaturePlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "percent.mito") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of mitochondrial genes or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.percent.mito") + 
           ggtitle( "Map of cells by filter status\n on %mitochondrial gene")
  )
}

# Plot the cells RNA counts in UMAP
FeaturePlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nCount_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with RNA counts")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by number of UMI count or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.nCount_RNA") + 
           ggtitle( "Map of cells by filter status \non nCount_RNA value")
  )
}

# Plot the cells genes counts in UMAP
FeaturePlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nFeature_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with \nnumber of detected genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by number of feature count or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.nFeature_RNA") + 
           ggtitle( "Map of cells by filter status \non nFeature_RNA value")
  )
}
 
# ..........................................................................................................
## @knitr dimreduc_ggplot_doublets
# ..........................................................................................................

# Plot the location and distribution of doublet cells in UMAP                                                     
DimPlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
         group.by = "scDblFinder") + 
  ggtitle( "Map of cells status on doublet")
