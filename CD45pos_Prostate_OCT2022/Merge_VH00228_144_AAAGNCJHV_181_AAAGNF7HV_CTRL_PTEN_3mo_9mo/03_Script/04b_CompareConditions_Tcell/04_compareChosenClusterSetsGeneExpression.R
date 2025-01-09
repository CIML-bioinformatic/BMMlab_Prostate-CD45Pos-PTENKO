# ##################################################
# This script analyses the expression of genes
# across clusters and samples
# ##################################################


# Compare expression of gene between sets of clusters 
# ######################################################

## @knitr compare_chosen_clusters_sets_expression

for( current_comparison in names( CHOSEN_CLUSTERS_COMPARISONS_LIST)){
  
  cluster_set_1 = CHOSEN_CLUSTERS_COMPARISONS_LIST[[ current_comparison]][[ "cluster_set_1"]]
  cluster_set_2 = CHOSEN_CLUSTERS_COMPARISONS_LIST[[ current_comparison]][[ "cluster_set_2"]]
  
  cat("\n \n")
  cat("### DEG analysis between clusters", paste( cluster_set_1, collapse = ";"), "and clusters", paste( cluster_set_2, collapse = ";"), "{.tabset .tabset-fade}")
  cat("\n \n")
  
  # Add a sub-tab with the DEG analysis results
  # .............................................
  
  cat("\n \n")
  cat("#### DEG analysis result\n")
  cat("\n \n")
  
  
  Idents( merge_seurat_object) = "seurat_clusters"
  nb_cell_in_cluster_set_1 = length( Idents( merge_seurat_object)[ which( Idents( merge_seurat_object) %in% cluster_set_1)])
  nb_cell_in_cluster_set_2 = length( Idents( merge_seurat_object)[ which( Idents( merge_seurat_object) %in% cluster_set_2)])
  
  all_cluster_comparison_markergenes_df = FindMarkers( merge_seurat_object, 
                                                           ident.1 = cluster_set_1, 
                                                           ident.2 = cluster_set_2,
                                                           test.use = DEG_TEST_USE)
  
  # Remove the possible infinite Log2FC values
  inf_log2fc_index = which( all_cluster_comparison_markergenes_df$avg_log2FC == "Inf" | all_cluster_comparison_markergenes_df$avg_log2FC == "-Inf")
  if( length( inf_log2fc_index) > 0){
    all_cluster_comparison_markergenes_df = all_cluster_comparison_markergenes_df[ -inf_log2fc_index,]
  }
  
  # Select genes with low enough p-value
  cluster_comparison_markergenes_df = all_cluster_comparison_markergenes_df[ which( all_cluster_comparison_markergenes_df$p_val_adj <= MARKER_GENES_ALPHA_THRESHOLD), ]
  
  if( nrow( cluster_comparison_markergenes_df) > 0){
    
    cluster_comparison_markergenes_set = row.names( cluster_comparison_markergenes_df)
    
    # Plot the volcano plot of the DEG
    print(
      EnhancedVolcano( toptable = all_cluster_comparison_markergenes_df,
                       lab = row.names( all_cluster_comparison_markergenes_df),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       selectLab = row.names( cluster_comparison_markergenes_df),
                       xlab = bquote(~Log[2]~ 'fold change'),
                       pCutoff = 0.05,
                       FCcutoff = 1.2,
                       legendPosition = "right",
                       legendLabSize = 8,
                       drawConnectors = TRUE,
                       widthConnectors = 0.75,
                       title = paste( "Clusters", paste( cluster_set_1, collapse = ";"), "vs", paste( cluster_set_2, collapse = ";")),
                       subtitle = paste( paste( cluster_set_1, collapse = ";"), "=", nb_cell_in_cluster_set_1, "cells",
                                  "/", paste( cluster_set_2, collapse = ";") ,"=", nb_cell_in_cluster_set_2, "cells")
      )
    )
    
    cluster_comparison_markergenes_df$gene = row.names( cluster_comparison_markergenes_df)
    cluster_comparison_markergenes_df$cluster_set_1 = paste( cluster_set_1, collapse = ";")
    cluster_comparison_markergenes_df$cluster_set_2 = paste( cluster_set_2, collapse = ";")
    row.names( cluster_comparison_markergenes_df) = NULL
    
    
    write.csv( cluster_comparison_markergenes_df,
               file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "clusterDEG_", paste( cluster_set_1, collapse = "-"), "vs", paste( cluster_set_2, collapse = "-"), "_", DEG_TEST_USE, ".csv")),
               row.names = cluster_comparison_markergenes_df$gene)
    
  }else{
    cluster_comparison_markergenes_set = NULL
    cat("<BR><b>No marker genes found for this comparison at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
  
  # Add a sub-tab with the enrichment of the positive markers in the GO BP terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n#### GO Biological Process enrichment in DEG")
  cat( "\n \n")
  
  # If some markers genes have been found, compute enrichment in GO BP
  if( ! is.null( cluster_comparison_markergenes_set)){
    
    # Compute enrichments
    positive_ego <- enrichGO(gene          = cluster_comparison_markergenes_set,
                             universe      = rownames( merge_seurat_object@assays$RNA),
                             OrgDb         = org.Mm.eg.db,
                             keyType       = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                             qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
    
    if( !is.null( positive_ego)){
      # Check for results
      selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
      
      # Plot the results if any, educing term by similarity of terms
      if( selected_term_number> 0){  
        
        print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of marker genes in BP GO terms"))
        positive_ego_df = data.frame( positive_ego)
        
        if( selected_term_number > 1){
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
      }
      else{
        cat("<BR>No result<BR>")
      }
    }else{
      cat("<BR><b>An issue occured during enrichment analysis</b>")
    }
  }else{
    cat("<BR><b>No marker genes found for this comparison at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
  
  
  
  # Add a sub-tab with the enrichment of the positive markers in the GO MF terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n#### GO Molecular Function enrichment in DEG")
  cat( "\n \n")
  
  # If some markers genes have been found, compute enrichment in GO MF
  if( ! is.null( cluster_comparison_markergenes_set)){
    
    # Compute enrichments
    positive_ego <- enrichGO(gene          = cluster_comparison_markergenes_set,
                             universe      = rownames( merge_seurat_object@assays$RNA),
                             OrgDb         = org.Mm.eg.db,
                             keyType       = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                             qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
    
    if( !is.null( positive_ego)){
      
      # Check for results
      selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
      
      # Plot the result if any, reducing by similarity of terms
      if( selected_term_number> 0){  
        
        print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of positive markers in MF GO terms"))
        positive_ego_df = data.frame( positive_ego)
        
        if( selected_term_number > 1){
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
      }
      else{
        cat("<BR>No result<BR>")
      }
    }else{
      cat("<BR><b>An issue occured during enrichment analysis</b>")
    }
  }else{
    cat("<BR><b>No marker genes found for cluster", clusterid, " at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
  
  
  cat( "\n \n")
  cat( "\n#### KEGG pathway enrichment in DEG")
  cat( "\n \n")
  
  # If some markers genes have been found, compute enrichment in GO MF
  if( ! is.null( cluster_comparison_markergenes_set)){
    
    # Get the entrezid of the Seurat object genes
    universe_entrez_id = AnnotationDbi::mapIds( org.Mm.eg.db, rownames( merge_seurat_object@assays$RNA), keytype="SYMBOL", column="ENTREZID")
    
    # Convert the marker gene symbols into entrezid
    markers_entrez_id = AnnotationDbi::mapIds( org.Mm.eg.db, cluster_comparison_markergenes_set, keytype="SYMBOL", column="ENTREZID")
    
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
    
    cat(" \n \n"); # Required for '.tabset'
  }else{
    cat("<BR><b>No marker genes found for this comparison at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
  
  cat(" \n \n"); # Required for '.tabset'
}




