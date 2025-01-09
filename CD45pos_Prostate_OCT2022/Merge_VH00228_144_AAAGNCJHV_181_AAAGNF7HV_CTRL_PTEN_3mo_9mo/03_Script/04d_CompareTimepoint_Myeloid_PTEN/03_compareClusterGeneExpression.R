# ##################################################
# This script analyses the expression of genes
# across clusters and samples
# ##################################################


# Compare expression of gene in cluster between samples
# ######################################################

## @knitr compare_cluster_expression

cat("\n \n")
cat("### DEG analysis between", DATE_NAME_9MO, "and", DATE_NAME_3MO, "in clusters {.tabset .tabset-fade}\n")

# Detect the DEG between samples in each cluster
# ...............................................
sample_cluster_markergenes_df = data.frame()

Idents( merge_seurat_object) = "seurat_clusters"
for( clusterid in levels( Idents( merge_seurat_object))){
  
  cat("\n \n")
  cat("#### Cluster", clusterid, " {.tabset .tabset-fade}\n");
  cat("\n \n")
  
  nb_cell_3mo_in_cluster = date_vs_cluster_df[ which( date_vs_cluster_df$cells_date == DATE_NAME_3MO & date_vs_cluster_df$cells_clusterid == clusterid), "Freq"]
  nb_cell_9mo_in_cluster = date_vs_cluster_df[ which( date_vs_cluster_df$cells_date == DATE_NAME_9MO & date_vs_cluster_df$cells_clusterid == clusterid), "Freq"]
  
  cat("<BR>")
  cat("<BR><b>Number of cell in cluster", clusterid, "from", DATE_NAME_3MO, ":", nb_cell_3mo_in_cluster, "</b>")
  cat("<BR><b>Number of cell in cluster", clusterid, "from", DATE_NAME_9MO, ":", nb_cell_9mo_in_cluster, "</b>")
  cat("<BR>")
  
  if( nb_cell_3mo_in_cluster < 10 || nb_cell_9mo_in_cluster < 10){
    cat("<BR><b>No statistics for this cluster because de nb of cells is too low (", nb_cell_3mo_in_cluster,"for 3mo cells,", nb_cell_9mo_in_cluster, "for 9mo cells)</b>")
    next
  }
  
  all_current_sample_cluster_markergenes_df = FindMarkers( merge_seurat_object, 
                                                           ident.1 = DATE_NAME_9MO, 
                                                           group.by = 'orig.date', 
                                                           subset.ident = clusterid,
                                                           test.use = DEG_TEST_USE)
  
  # Remove the possible infinite Log2FC values
  inf_log2fc_index = which( all_current_sample_cluster_markergenes_df$avg_log2FC == "Inf" | all_current_sample_cluster_markergenes_df$avg_log2FC == "-Inf")
  if( length( inf_log2fc_index) > 0){
    all_current_sample_cluster_markergenes_df = all_current_sample_cluster_markergenes_df[ -inf_log2fc_index,]
  }
  
  
  all_current_sample_cluster_markergenes_df$diff.pct = all_current_sample_cluster_markergenes_df$pct.1 - all_current_sample_cluster_markergenes_df$pct.2
  all_current_sample_cluster_markergenes_df = all_current_sample_cluster_markergenes_df[ order( all_current_sample_cluster_markergenes_df$diff.pct, decreasing = TRUE), ]
  
  # Select genes with low enough p-value
  current_sample_cluster_markergenes_df = all_current_sample_cluster_markergenes_df[ which( all_current_sample_cluster_markergenes_df$p_val_adj <= MARKER_GENES_ALPHA_THRESHOLD), ]
  
  if( nrow( current_sample_cluster_markergenes_df) > 0){
    
    # Accumulate information along clusters
    current_sample_cluster_markergenes_df = current_sample_cluster_markergenes_df[ order( current_sample_cluster_markergenes_df$diff.pct, decreasing = TRUE), ]
    current_sample_cluster_markergenes_df$gene = row.names( current_sample_cluster_markergenes_df)
    current_sample_cluster_markergenes_df$cluster = clusterid
    row.names( current_sample_cluster_markergenes_df) = NULL
    sample_cluster_markergenes_df = rbind( sample_cluster_markergenes_df, current_sample_cluster_markergenes_df)
    
    
    # Plot the Volcano plot of the results
    # .....................................
    
    cat("\n \n")
    cat("##### Volcano Plot\n")
    cat("\n \n")
    
    # Plot the volcano plot of the DEG
    print(
      EnhancedVolcano( toptable = all_current_sample_cluster_markergenes_df,
                       lab = row.names( all_current_sample_cluster_markergenes_df),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       selectLab = row.names( current_sample_cluster_markergenes_df),
                       xlab = bquote(~Log[2]~ 'fold change'),
                       pCutoff = MARKER_GENES_ALPHA_THRESHOLD,
                       FCcutoff = MARKER_GENES_LOG2FC_THRESHOLD,
                       legendPosition = "right",
                       legendLabSize = 3,
                       drawConnectors = TRUE,
                       widthConnectors = 0.75,
                       title = paste( "DEG beetween", DATE_NAME_9MO, "and", DATE_NAME_3MO),
                       subtitle = paste( "Cluster", clusterid, ": 9mo=", nb_cell_9mo_in_cluster, "cells / 3mo=", nb_cell_3mo_in_cluster, "cells")
      )
    )
    
    
    # Plot the Feature plots of the results
    # .....................................
    
    # cat(" \n \n")
    # cat("##### Marker genes expression projection")
    # cat(" \n \n"); # Required for '.tabset'
    # 
    # # Plot expression values of marker genes on dimreduc figures for the current cluster
    # 
    # # Highlight cells of current cluster on a dimreduc plot
    # highlightClusterPlot(clusterid, seuratObject = merge_seurat_object, reduction = ifelse( exists("useReduction"), useReduction, "umap"));
    # 
    # # Plots expression on projected cells
    # invisible( lapply( current_sample_cluster_markergenes_df$gene[ 1:9], function(featureName)
    # {
    #   print( FeaturePlot( merge_seurat_object, 
    #                       features = featureName, 
    #                       reduction = ifelse( exists("useReduction"), useReduction, "umap"), 
    #                       order = TRUE,
    #                       split.by = "orig.date") +
    #            theme( axis.title.x = element_blank(),
    #                   axis.title.y = element_blank(),
    #                   legend.position = "none"));
    # }));
    
    cat(" \n \n"); # Required for '.tabset'
    
    # Compute the GO enrichment for each list of marker genes of the current cluster

    # Add a sub-tab with the enrichment of the positive markers in the GO BP terms
    # ...............................................................................
    cat( "\n \n")
    cat( "#####  Markers genes GO Biological Process (positive markers)")
    cat( "\n \n")
    
    positive_current_sample_cluster_markergenes_df = current_sample_cluster_markergenes_df[ which( current_sample_cluster_markergenes_df$avg_log2FC > MARKER_GENES_LOG2FC_THRESHOLD), ]
    
    if( nrow( positive_current_sample_cluster_markergenes_df) > 30){
      
      # Compute enrichments
      positive_ego <- enrichGO(gene          = positive_current_sample_cluster_markergenes_df$gene[ 1:min( 200, nrow( positive_current_sample_cluster_markergenes_df))],
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
        
        # Save result to file
        positive_ego_df = data.frame( positive_ego)
        positive_ego_df$cluster = clusterid
        write.table( positive_ego_df, col.names = TRUE, row.names = FALSE, sep=";",
                     file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "GOBPEnrichment_PositiveMarkers_Cluster", clusterid, ".csv"))
        )
        
        # Plot the results
        print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of marker genes in BP GO terms (positive)"))
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
    }
    else{
      cat("<BR><b>Not enough positive marker genes found for cluster", clusterid, "</b>")
    }
    
    cat(" \n \n"); # Required for '.tabset'
    
    # Add a sub-tab with the enrichment of the negative markers in the GO BP terms
    # ...............................................................................
    cat( "\n \n")
    cat( "#####  Markers genes GO Biological Process (negative markers)")
    cat( "\n \n")
    
    negative_current_sample_cluster_markergenes_df = current_sample_cluster_markergenes_df[ which( current_sample_cluster_markergenes_df$avg_log2FC < MARKER_GENES_LOG2FC_THRESHOLD), ]
    
    if( nrow( negative_current_sample_cluster_markergenes_df) > 30){
      
      # Compute enrichments
      negative_ego <- enrichGO(gene          = negative_current_sample_cluster_markergenes_df$gene[ 1:min( 200, nrow( negative_current_sample_cluster_markergenes_df))],
                               universe      = rownames( merge_seurat_object@assays$RNA),
                               OrgDb         = org.Mm.eg.db,
                               keyType       = "SYMBOL",
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                               qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
      
      # Check for results
      selected_term_number = length( which( negative_ego@result$p.adjust <= negative_ego@pvalueCutoff & negative_ego@result$qvalue <= negative_ego@qvalueCutoff))
      
      # Plot the results if any, educing term by similarity of terms
      if( selected_term_number> 0){  
        
        # Save result to file
        negative_ego_df = data.frame( negative_ego)
        negative_ego_df$cluster = clusterid
        write.table( negative_ego_df, col.names = TRUE, row.names = FALSE, sep=";",
                     file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "GOBPEnrichment_NegativeMarkers_Cluster", clusterid, ".csv"))
        )
        
        # Plot the results
        print( dotplot( negative_ego, showCategory=10) + ggtitle("Enrichment of marker genes in BP GO terms (negative)"))
        negative_ego_df = data.frame( negative_ego)
        
        negative_ego_simMatrix <- calculateSimMatrix( negative_ego_df$ID,
                                                      orgdb="org.Mm.eg.db",
                                                      ont="BP",
                                                      method="Rel")
        
        negative_ego_scores <- setNames(-log10( negative_ego_df$qvalue), negative_ego_df$ID)
        negative_ego_reducedTerms <- reduceSimMatrix(negative_ego_simMatrix,
                                                     negative_ego_scores,
                                                     threshold=0.7,
                                                     orgdb="org.Mm.eg.db")
        
        print( scatterPlot( negative_ego_simMatrix, negative_ego_reducedTerms))
        
        treemapPlot( negative_ego_reducedTerms)
      }
      else{
        cat("<BR>No result<BR>")
      }
    }
    else{
      cat("<BR><b>Not enough negative marker genes found for cluster", clusterid, "</b>")
    }
    
    cat(" \n \n"); # Required for '.tabset'
    
  }else{
    cat("<BR><b>No marker genes found for cluster", clusterid, " at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
}

cat(" \n \n"); # Required for '.tabset'

all_sample_cluster_markergenes = unique( sample_cluster_markergenes_df$gene)

write.csv( sample_cluster_markergenes_df,
           file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "clusterDEG_", DATE_NAME_9MO, "vs", DATE_NAME_3MO, "_", DEG_TEST_USE, ".csv")),
           row.names = sample_cluster_markergenes_df$gene)


