# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

# MARKER GENES
##############

# .............................................................
## @knitr heterogeneity_markerGenes
# .............................................................

# Identify marker genes
Idents( merge_seurat_object) = "seurat_clusters"
markers = FindAllMarkers( object          = merge_seurat_object,
                          test.use        = FINDMARKERS_METHOD,
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          return.thresh   = FINDMARKERS_PVAL_THR,
                          random.seed     = SEED,
                          verbose         = .VERBOSE);

# Save markers list as 'tsv' table
write.table( markers,
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "MarkerGenes.tsv")),
             quote = FALSE,
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t")

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers_table = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP_TABLE)) x else head( x, n = FINDMARKERS_SHOWTOP_TABLE));
});

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers_heatmap = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP_HEATMAP)) x else head( x, n = FINDMARKERS_SHOWTOP_HEATMAP));
});

# Merge marker genes in a single data.frame and render it as datatable
topMarkers_tableDF = do.call( rbind, topMarkers_table);
topMarkers_heatmapDF = do.call( rbind, topMarkers_heatmap);

# Select and order columns to be shown in datatable
topMarkers_tableDT = topMarkers_tableDF[c("gene", "cluster", "avg_log2FC", "p_val_adj")]



# .............................................................
## @knitr heterogeneity_markerGenes_table
# .............................................................

cat(" \n \n")
cat("### Marker genes datatable {.tabset .tabset-fade}")
cat(" \n \n") # Required for '.tabset'

# Create datatable
datatable( topMarkers_tableDT,
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Avg. LogFC", "Adj. Pvalue"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP_TABLE), "All", paste("Top", FINDMARKERS_SHOWTOP_TABLE)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'))%>%
  # Add bar relative to logFC
  formatStyle( columns = "avg_log2FC",
               background = styleColorBar( data = range( topMarkers_tableDT[["avg_log2FC"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  # Add color from cluster
  formatStyle( columns = "cluster",
               backgroundColor = styleEqual( names(clustersColor),
                                             scales::alpha(clustersColor, 0.3)));



# .............................................................
## @knitr heterogeneity_markerGenes_heatmap_mean
# .............................................................

cat(" \n \n")
cat("### Marker genes heatmaps/dotplot {.tabset .tabset-fade}")
cat(" \n \n") # Required for '.tabset'

# Get the matrix of expression and associated clusters from Seurat object
# ........................................................................
expMat = as.matrix( LayerData( merge_seurat_object, assay = "RNA", layer = "data"));
Idents( merge_seurat_object) = "seurat_clusters"
clusterID = Idents( merge_seurat_object);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkers_heatmapDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Compute the mean of top markers in clusters and produce matrix with the result
# ........................................................................
all_mean_expression = vector()
cluster_set = unique( topMarkers_heatmapDF$cluster)
clusters_marker_genes = topMarkers_heatmapDF[ which( topMarkers_heatmapDF$cluster %in% cluster_set), "gene"]
for( cluster_id in cluster_set){
  clusters_cells = names( clusterID)[ which( clusterID == cluster_id)]
  all_mean_expression = c( all_mean_expression, BiocGenerics::rowMeans( expMat[ clusters_marker_genes, clusters_cells]))
}
meanExpMat = t( matrix( all_mean_expression, byrow = TRUE, ncol = FINDMARKERS_SHOWTOP_HEATMAP*length( cluster_set)))
colnames( meanExpMat) = cluster_set
rownames( meanExpMat) = paste0( topMarkers_heatmapDF[ which( topMarkers_heatmapDF$cluster %in% cluster_set), "cluster"], '.', clusters_marker_genes)

# Plot a heatmap of the matrix of mean expression of top marker genes in clusters
# ............................................................................

cat("\n \n")
pheatmap( expMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(80),
          breaks = seq( -4.2, 4.2, 0.1),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkers_heatmapDF[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = clusterID, stringsAsFactors = FALSE, row.names = colnames( expMat)),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = FALSE,
          fontsize_row = 5,
          main = "Scaled normalized expression by cell\nof top markers");

cat("\n \n")
Idents( merge_seurat_object) = "seurat_clusters"
DotPlot( merge_seurat_object, 
         assay = "RNA",
         features = unique( clusters_marker_genes),
         cols = "RdBu") +
  theme( axis.text.x = element_text( angle = 45, hjust = 1)) +
  ggtitle( "Scaled normalized mean expression of topmarker genes in clusters")

cat("\n \n")
pheatmap( meanExpMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkers_heatmapDF[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = cluster_set, stringsAsFactors = FALSE, row.names = cluster_set),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = TRUE,
          fontsize_row = 5,
          main = "Scaled mean normalized expression\nby cluster of top markers");


# .............................................................
## @knitr heterogeneity_markerGenes_expression_projection
# .............................................................

cat(" \n \n")
cat("### Marker genes expression projection {.tabset .tabset-fade}")
cat(" \n \n"); # Required for '.tabset'

# Plot expression values of marker genes on dimreduc figures for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers_heatmap), function(clusterName)
{
  cat("#### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");

  # Highlight cells of current cluster on a dimreduc plot
  highlightClusterPlot(clusterName, seuratObject = merge_seurat_object, reduction = ifelse( exists("useReduction"), useReduction, "umap"));

  # Plots expression on projected cells
  invisible( lapply( topMarkers_heatmap[[clusterName]][["gene"]], function(featureName)
    {
      print( FeaturePlot( merge_seurat_object, features = featureName, reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
               theme( axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      legend.position = "none"));
    }));

  cat(" \n \n"); # Required for '.tabset'
}));



# .............................................................
## @knitr heterogeneity_markerGenes_expression_violin
# .............................................................

cat("### Marker genes expression violin plot {.tabset .tabset-fade}")
cat(" \n \n"); # Required for '.tabset'

# Plot expression values of marker genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers_heatmap), function(clusterName)
{
  cat("#### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");

  # Remind cluster name in an empty figure to keep consistent alignment of panels between tabs
  plot( c( 0, 1), c( 0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
  text( x = 0.5, y = 0.5, paste( "Cluster", clusterName), cex = 2, col = clustersColor[clusterName]);

  # Violinplot for expression value of marker genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( topMarkers_heatmap[[clusterName]][["gene"]], violinFeatureByCluster, seuratObject = merge_seurat_object, clustersColor = clustersColor));

  cat(" \n \n"); # Required for '.tabset'
}))


# .............................................................
## @knitr heterogeneity_markerGenes_functional_enrichment
# .............................................................

cat(" \n \n")
cat("### Marker genes functional enrichment {.tabset .tabset-fade}")
cat(" \n \n") # Required for '.tabset'

# Compute the GO enrichment for each list of marker genes of clusters
invisible( lapply( names( topMarkers_heatmap), function( clusterName)
{
  cat( "\n \n")
  cat("#### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[ clusterName], "; padding:0px 2px'>", clusterName, "</span>  {.tabset .tabset-fade}");
  cat( "\n \n")
  
  # Add a sub-tab with the enrichment of the positive markers in the GO BP terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### Markers genes GO Biological Process")
  cat( "\n \n")
  
  # Compute enrichments
  positive_ego <- enrichGO(gene          = markers[ which( markers$cluster == clusterName), "gene"],
                           universe      = rownames( merge_seurat_object@assays$RNA),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                           qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
  
  # Check for results
  selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
  
  # Plot the results if any, educing term by similarity of terms
  if( selected_term_number> 0){  
    
    print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of marker genes in BP GO terms"))
    positive_ego_df = data.frame( positive_ego)
    
    positive_ego_simMatrix <- calculateSimMatrix( positive_ego_df$ID,
                                                  orgdb="org.Mm.eg.db",
                                                  ont="BP",
                                                  method="Rel")
    
    positive_ego_scores <- setNames(-log10( positive_ego_df$qvalue), positive_ego_df$ID)
    positive_ego_reducedTerms <- reduceSimMatrix(positive_ego_simMatrix,
                                                 positive_ego_scores,
                                                 threshold=0.7,
                                                 orgdb="org.Mm.eg.db")
    
    print( scatterPlot( positive_ego_simMatrix, positive_ego_reducedTerms))
    
    treemapPlot( positive_ego_reducedTerms)
  }
  else{
    cat("<BR>No result<BR>")
  }
  
  
    # Add a sub-tab with the enrichment of the positive markers in the GO MF terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### Markers genes GO Molecular Function")
  cat( "\n \n")
  
  # Compute enrichments
  positive_ego <- enrichGO(gene          = markers[ which( markers$cluster == clusterName), "gene"],
                           universe      = rownames( merge_seurat_object@assays$RNA),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                           qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
  
  # Check for results
  selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
  
  # Plot the result if any, reducing by similarity of terms
  if( selected_term_number> 0){  
    
    print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of positive markers in MF GO terms"))
    positive_ego_df = data.frame( positive_ego)
    
    positive_ego_simMatrix <- calculateSimMatrix( positive_ego_df$ID,
                                                  orgdb="org.Mm.eg.db",
                                                  ont="MF",
                                                  method="Rel")
    
    positive_ego_scores <- setNames(-log10( positive_ego_df$qvalue), positive_ego_df$ID)
    positive_ego_reducedTerms <- reduceSimMatrix(positive_ego_simMatrix,
                                                 positive_ego_scores,
                                                 threshold=0.7,
                                                 orgdb="org.Mm.eg.db")
    
    print( scatterPlot( positive_ego_simMatrix, positive_ego_reducedTerms))
    
    treemapPlot( positive_ego_reducedTerms)
  }
  else{
    cat("<BR>No result<BR>")
  }

  cat(" \n \n"); # Required for '.tabset'
  cat("\n##### Marker genes KEGG pathways")
  cat(" \n \n"); # Required for '.tabset'
  
  # Get the entrez id of the Seurat object genes
  universe_entrez_id = AnnotationDbi::mapIds( org.Mm.eg.db, rownames( merge_seurat_object@assays$RNA), keytype="SYMBOL", column="ENTREZID")

  # Convert the marker gene symbols into entrezid
  markers_entrez_id = AnnotationDbi::mapIds( org.Mm.eg.db, markers[ which( markers$cluster == clusterName), "gene"], keytype="SYMBOL", column="ENTREZID")
  
  # Compute the enrichment
  positive_kegg <- enrichKEGG( gene         = markers_entrez_id,
                               organism      = "mmu",
                               universe      = universe_entrez_id,
                               keyType       = "kegg",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                               qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
  
  if( !is.null( positive_kegg)){
    positive_kegg_result_df = positive_kegg@result
  
    # Plot the results for all KEGG Pathways
    
    cat("<BR><b>Enrichment considering all KEGG pathways</b>")
    if( sum( positive_kegg_result_df$p.adjust <= ENRICHMENT_GO_PVALUECUTOFF) > 0){  
      print( dotplot( positive_kegg))
    }else{
      cat("<BR>No enrichment found")
    }
  
    # Plot the results for the chosen KEGG Pathways (recomputing the mutitest adjustment with limited list)
    
    cat("\n<b>Enrichment considering only the chosen KEGG pathways</b>")
    selected_positive_kegg_result_df = positive_kegg_result_df[ which( positive_kegg_result_df$ID %in% KEGG_PATHWAY_DF$Pathway_id), ]
    selected_positive_kegg_result_df$p.adjust = p.adjust( selected_positive_kegg_result_df$pvalue)
    selected_positive_kegg_result_df$qvalue = p.adjust( selected_positive_kegg_result_df$pvalue)
    if( sum( selected_positive_kegg_result_df$p.adjust <= ENRICHMENT_GO_PVALUECUTOFF) > 0){
      positive_kegg@result = selected_positive_kegg_result_df
      print( dotplot( positive_kegg))
    }else{
      cat("<BR>No enrichment found")
    }
  }else{
    cat("<BR>Enrichment analysis got an issue.")
    cat("<BR>markers_entrez_id = ", paste( markers_entrez_id, collapse = ";"))
  }
  
}))
