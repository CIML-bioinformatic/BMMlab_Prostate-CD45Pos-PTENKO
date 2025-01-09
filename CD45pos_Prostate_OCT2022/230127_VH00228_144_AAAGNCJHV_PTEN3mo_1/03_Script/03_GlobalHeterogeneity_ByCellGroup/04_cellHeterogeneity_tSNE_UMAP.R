# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# DIMENSIONAL REDUCTION (TSNE/UMAP)
###################################

# ..........................................................................................................
## @knitr heterogeneity_dimReduc
# ..........................................................................................................

nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
if(DIMREDUC_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_dimreduc = nbPC
}

sc10x = RunUMAP( sc10x, dims = 1:nbPC_dimreduc);
sc10x = RunTSNE( sc10x, dims = 1:nbPC_dimreduc);

# Save resulting coordinates for all cells as 'tsv' files
write.table( Embeddings(sc10x, reduction = "umap"), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_umap.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

write.table( Embeddings(sc10x, reduction = "tsne"), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_tsne.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");


# ..........................................................................................................
## @knitr dimreduc_ggplot_covariables
# ..........................................................................................................

# Plot the cells sample of origin in UMAP
DimPlot( sc10x, reduction = useReduction, group.by = "original_seurat_clusters") + 
  ggtitle( "Map of cells by cluster of origin")


# Plot the cells percentage of ribosomal genes in UMAP
FeaturePlot( sc10x, features = "percent.ribo") +
              scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
              ggtitle( "Map of cells with level of\n percentage of ribosomal genes")

# Plot the cells percentage of mitochondrial genes in UMAP  
FeaturePlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "percent.mito") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")

# Plot the cells RNA counts in UMAP
FeaturePlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nCount_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with RNA counts")

# Plot the cells genes counts in UMAP
FeaturePlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nFeature_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with \nnumber of detected genes")


