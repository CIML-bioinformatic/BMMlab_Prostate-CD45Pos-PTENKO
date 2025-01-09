# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# IDENTIFY CLUSTERS
###################

# ..........................................................................................................
## @knitr heterogeneity_identifyClusters_makeClusters
# ..........................................................................................................

# Identify clusters of cells by graph approach
nbPC_findclusters=FINDCLUSTERS_USE_PCA_NBDIMS
if(FINDCLUSTERS_USE_PCA_NBDIMS>nbPC){
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'findclusters' (", FINDCLUSTERS_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_findclusters = nbPC
}  

merge_seurat_object <- FindNeighbors(object    = merge_seurat_object,
                                     k.param   = FINDNEIGHBORS_K, 
                                     reduction = "pca",
                                     dims      = 1:nbPC_findclusters,
                                     verbose   = .VERBOSE);


# Compute the clusters along several resolutions
for( resolution in FINDCLUSTERS_RESOLUTION_RANGE){
  
  merge_seurat_object = FindClusters( object             = merge_seurat_object,
                                      resolution         = resolution,
                                      algorithm          = FINDCLUSTERS_ALGORITHM,
                                      temp.file.location = "/tmp/",
                                      verbose            = FALSE);
}

# Display clusters along various resolution with clustree
# ....................................................................

cat(" \n \n")
cat("### Clustree map")
cat(" \n \n")

clustree( merge_seurat_object, prefix = "RNA_snn_res.", layout = "sugiyama")
clustree( merge_seurat_object, prefix = "RNA_snn_res.", layout = "sugiyama", node_colour = "sc3_stability")

# Display the various clustering levels in UMAP
# ....................................................................

invisible( lapply( FINDCLUSTERS_RESOLUTION_RANGE, function( resolution)
{
  
  cat(" \n \n")
  cat("### Resolution", resolution)
  cat(" \n \n")
  
  Idents( merge_seurat_object) = paste0( "RNA_snn_res.", resolution)
  cluster_number = length( levels( Idents( merge_seurat_object)))
  print(
    tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
      DimPlot( merge_seurat_object, reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
               order = TRUE, 
               label = TRUE, 
               label.size = 6)  +
        ggtitle( paste0( "Resolution ", resolution, " (", cluster_number, " clusters)")) + 
        theme( axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "none")),
      error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  
  cat(" \n \n"); # Required for '.tabset'
}));


# Choose the resolution and display the cluster size at the chosen resolution
# .............................................................................

cat(" \n \n")
cat("### Cluster size at chosen resolution (", FINDCLUSTERS_RESOLUTION, ")")
cat(" \n \n")

# Choose the cluster resolution to be used among all those computed
if( ! FINDCLUSTERS_RESOLUTION %in% FINDCLUSTERS_RESOLUTION_RANGE){
  stop( paste( "The cluster resolution", FINDCLUSTERS_RESOLUTION, "is not is the studied range of resolution:", paste( FINDCLUSTERS_RESOLUTION_RANGE, collapse = ";")))
}
Idents( merge_seurat_object) = paste0( "RNA_snn_res.", FINDCLUSTERS_RESOLUTION)
merge_seurat_object = AddMetaData( merge_seurat_object, metadata = Idents( merge_seurat_object), col.name = "seurat_clusters")

# Rename the clusters if required
if( exists( "RENAME_CLUSTERS") && SAMPLE_NAME %in% names( RENAME_CLUSTERS)){
  Idents( merge_seurat_object) = "seurat_clusters"
  merge_seurat_object = AddMetaData( merge_seurat_object, metadata = Idents( merge_seurat_object), col.name = "original_seurat_clusters")
  new_names = Idents( merge_seurat_object)
  for( cluster_id in levels( Idents( merge_seurat_object))){
    new_names[ which( Idents( merge_seurat_object) == cluster_id)] = RENAME_CLUSTERS[[ SAMPLE_NAME]][[ cluster_id]]
  }
  merge_seurat_object = AddMetaData( merge_seurat_object, metadata = new_names, col.name = "seurat_clusters")
  
}

# Show number of cells in each cluster
clustersCount = as.data.frame( table( Cluster = merge_seurat_object[[ "seurat_clusters" ]]), responseName = "CellCount");

# Define a set of colors for clusters (based on ggplot default)
Idents( merge_seurat_object) = "seurat_clusters"
clustersColor = scales::hue_pal()( base::nlevels( Idents( merge_seurat_object)))
names( clustersColor) = levels( Idents( merge_seurat_object))


# Save cells cluster identity as determined with 'FindClusters'
write.table( data.frame( merge_seurat_object[[ "numID"]], identity = Idents(merge_seurat_object)), 
             file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "cellsClusterIdentity.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t")

# Also save cluster color attribution for reference
# Save cells cluster identity as determined with 'FindClusters'
write.table( clustersColor, 
             file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "clustersColor.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = FALSE,
             sep="\t")


# Create datatable
datatable( clustersCount,
           class = "compact",
           rownames = FALSE,
           colnames = c("Cluster", "Nb. Cells"),
           options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          columnDefs = list( # Center all columns
                                            list( targets = 0:(ncol(clustersCount)-1),
                                            className = 'dt-center')),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          paging = FALSE, # Disable pagination (show all)
                          processing = TRUE,
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%
  # Add color from cluster
  formatStyle( columns = "seurat_clusters",
               color = styleEqual( names(clustersColor), clustersColor),
               fontWeight = 'bold')


# Compute a matrix of average expression value by cluster (for each gene)
geneExpByCluster = do.call( rbind, 
                            apply( as.matrix( LayerData( merge_seurat_object)), # expression values
                                   1,                                # by rows
                                   tapply,                           # apply by group
                                   INDEX = Idents( merge_seurat_object),           # clusters IDs
                                   mean,                             # summary function
                                   simplify = FALSE));               # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByCluster, 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "normExpressionByCluster.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");



# ..........................................................................................................
## @knitr heterogeneity_identifyClusters_splitStats
# ..........................................................................................................

# Show UMIs Genes Mitochondrial and Ribosomal content split by cluster

# Gather data to be visualized together (cell name + numID + metrics)
Idents( merge_seurat_object) = "seurat_clusters"
cellsData = cbind( "Cell" = colnames( merge_seurat_object), # cell names from rownames conflict with search field in datatable
                   merge_seurat_object[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( merge_seurat_object$percent.mito)) as.numeric(format(merge_seurat_object[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( merge_seurat_object$percent.ribo)) as.numeric(format(merge_seurat_object[["percent.ribo", drop = TRUE]], digits = 5)) else NULL,
                   "Batch" = merge_seurat_object[["orig.ident"]],
                   "Cluster" = Idents( merge_seurat_object));


# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                     "Cell ID: ",
                                     "# UMIs: ",
                                     "# Genes: ",
                                     if(length( merge_seurat_object$percent.mito)) "% Mito: ",
                                     if(length( merge_seurat_object$percent.ribo)) "% Ribo: "),
                                  cellsData[-ncol( cellsData)],
                                  sep = ""),
                    sep = "\n"));

# Define size for panels (or assembled figure when using subplot)
panelWidth = 90 * nlevels(cellsData[["Cluster"]]);
panelHeight = 800;

# Generate plotly violin/jitter panels for #umis, #genes, and %mitochondrial stats
lypanel_umis  = plotViolinJitter( cellsData,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nCount_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# UMIs",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = panelWidth,
                                  panelHeight = panelHeight);

lypanel_genes = plotViolinJitter( cellsData,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nFeature_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# Genes",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = panelWidth,
                                  panelHeight = panelHeight);

lypanel_mitos = if(length( merge_seurat_object$percent.mito)) plotViolinJitter( cellsData,
                                                          xAxisFormula = ~as.numeric( Cluster),
                                                          yAxisFormula = ~percent.mito,
                                                          colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                          yAxisTitle = "% Mito",
                                                          hoverText = hoverText,
                                                          traceName = ~paste( "Cluster", Cluster),
                                                          xTicklabels = levels(cellsData[["Cluster"]]),
                                                          panelWidth = panelWidth,
                                                          panelHeight = panelHeight) else NULL;

lypanel_ribos = if(length( merge_seurat_object$percent.ribo)) plotViolinJitter( cellsData,
                                                          xAxisFormula = ~as.numeric( Cluster),
                                                          yAxisFormula = ~percent.ribo,
                                                          colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                          yAxisTitle = "% Ribo",
                                                          hoverText = hoverText,
                                                          traceName = ~paste( "Cluster", Cluster),
                                                          xTicklabels = levels(cellsData[["Cluster"]]),
                                                          panelWidth = panelWidth,
                                                          panelHeight = panelHeight) else NULL;

# Set panels as a list, define plotly config, and remove legend
panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
panelsList = lapply( panelsList, config, displaylogo = FALSE,
                     toImageButtonOptions = list( format='svg'),
                     modeBarButtons = list( list('toImage'),
                                            list( 'zoom2d', 'pan2d', 'resetScale2d')));

# Group plotly violin/jitter panels so we can synchronise axes and use highlight on the full set
# plotPanels = layout( subplot( panelsList,
#                               nrows = 4,
#                               shareX = TRUE,
#                               titleY = TRUE),
#                      xaxis = list(title = "Seurat Cluster",
#                                   showgrid = TRUE,
#                                   tickvals = seq(nlevels(cellsData[["Cluster"]])),
#                                   ticktext = levels(cellsData[["Cluster"]])),
#                      showlegend = FALSE, # Remove eventual legends (does not mix well with subplot and highlight)
#                      autosize = TRUE);

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
# div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
#      div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
#           div(plotPanels, style = paste("flex : 0 0 auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
#)


## @knitr heterogeneity_dimReduc_with_clusters

# Plot the map with clusters with ggplot
print(
DimPlot( merge_seurat_object, 
         reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
         group.by = "seurat_clusters",
         label = TRUE) + 
  ggtitle( "Map of cells with clusters")
)

# Plot the map with clusters with ggplot vs sample of origin
print(
  DimPlot( merge_seurat_object, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "seurat_clusters", split.by = "orig.ident", ncol = 2,
           label = TRUE) + 
    ggtitle( "Map of cells with clusters")
)

origin_vs_cluster_table = table( merge_seurat_object$orig.ident, merge_seurat_object$seurat_clusters)
origin_vs_cluster_table %>% kable(caption = "Cell distribution of Sample versus Clusters") %>% kable_styling()

originsample_vs_cluster_chisq_test = chisq.test( origin_vs_cluster_table)
corrplot::corrplot( originsample_vs_cluster_chisq_test$residuals, is.cor = FALSE)

# Plot the map with clusters with ggplot vs condition of origin
print(
  DimPlot( merge_seurat_object, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "seurat_clusters", split.by = "orig.condition", ncol = 2,
           label = TRUE) + 
    ggtitle( "Map of cells with clusters")
)

condition_vs_cluster_table = table( merge_seurat_object$orig.condition, merge_seurat_object$seurat_clusters)
condition_vs_cluster_table %>% kable( caption = "Cell distribution of Condition versus Clusters") %>% kable_styling()

origincondition_vs_cluster_chisq_test = chisq.test( condition_vs_cluster_table)
corrplot::corrplot( origincondition_vs_cluster_chisq_test$residuals, is.cor = FALSE)

