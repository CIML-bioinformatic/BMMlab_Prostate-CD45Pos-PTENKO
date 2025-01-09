# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# DIMENSIONAL REDUCTION (TSNE/UMAP)
###################################

# ..........................................................................................................
## @knitr heterogeneity_dimReduc
# ..........................................................................................................
useReduction = "umap"

nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
if(DIMREDUC_USE_PCA_NBDIMS>nbPC){
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_dimreduc = nbPC
}

# Using scDEED to obtain the best UMAP parameters
if( OPTIMIZE_UMAP_PARAMETERS == TRUE){
  cat("<BR><br>scDEED analysis</b><BR>")
  umap_analysis <- scDEED( merge_seurat_object@assays$RNA@counts,
                           num_pc = nbPC_dimreduc,
                           n_neighbors = c(seq(from=8,to=12,by=2)),
                           min.dist = seq(0.05,0.15, by = 0.05),
                           use_method = useReduction,
                           visualization = TRUE)
  
  best_UMAP_n.neighbors = umap_analysis$`best pair of n.neighbors and min.dist`[ 1, "n.neighbers"]
  best_UMAP_min.dist = umap_analysis$`best pair of n.neighbors and min.dist`[ 1, "min.dist"]

  cat("<BR><BR>")
}else{
  cat("<BR><br>Bypassing scDEED analysis</b><BR>")
  best_UMAP_n.neighbors = BEST_UMAP_N_NEIGHBORS 
  best_UMAP_min.dist = BEST_UMAP_MIN_DIST
}

cat("<BR>Best Neighbors number for UMAP =", best_UMAP_n.neighbors)
cat("<BR>Best Min distance number for UMAP =", best_UMAP_min.dist)

# Run the UMAP with the optimized parameters
merge_seurat_object = RunUMAP( merge_seurat_object, dims = 1:nbPC_dimreduc, n.neighbors = best_UMAP_n.neighbors, min.dist = best_UMAP_min.dist)

# Save resulting coordinates for all cells as 'tsv' files
write.table( Embeddings(merge_seurat_object, reduction = "umap"), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_umap.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# ..........................................................................................................
## @knitr dimreduc_ggplot_covariables
# ..........................................................................................................

# Plot the cells sample of origin in UMAP
DimPlot( merge_seurat_object, reduction = useReduction, group.by = "orig.ident") + 
  ggtitle( "Map of cells by sample of origin")


# if QC EXPLORATION MODE is active, plot if the cells are outliers or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
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
DimPlot( merge_seurat_object, reduction = useReduction, group.by = "scDblFinder") + 
  ggtitle( "Map of cells status on doublet")
