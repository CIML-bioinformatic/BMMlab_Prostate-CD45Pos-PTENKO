# ##################################################
# This script analyses the expression of genes
# across clusters and samples
# ##################################################


# Compare expression of gene in cluster between samples
# ######################################################

## @knitr compare_cluster_expression

cat("\n \n")
cat("### DEG analysis between", CONDITION_NAME_PTEN, "and", CONDITION_NAME_CTRL, "in clusters {.tabset .tabset-fade}\n")

# Detect the DEG between samples in each cluster
# ...............................................
sample_cluster_markergenes_df = data.frame()

Idents( merge_seurat_object) = "seurat_clusters"
for( clusterid in levels( Idents( merge_seurat_object))){
  
  cat("\n \n")
  cat("#### Cluster", clusterid, " {.tabset .tabset-fade}")
  cat("\n \n")
  
  cat("\n \n")
  cat("##### DEG analysis result\n")
  cat("\n \n")
  
  nb_cell_PTEN_in_cluster = condition_vs_cluster_df[ which( condition_vs_cluster_df$cells_condition == CONDITION_NAME_PTEN & condition_vs_cluster_df$cells_clusterid == clusterid), "Freq"]
  nb_cell_CTRL_in_cluster = condition_vs_cluster_df[ which( condition_vs_cluster_df$cells_condition == CONDITION_NAME_CTRL & condition_vs_cluster_df$cells_clusterid == clusterid), "Freq"]
  
  all_current_sample_cluster_markergenes_df = FindMarkers( merge_seurat_object, 
                                                           ident.1 = CONDITION_NAME_PTEN, 
                                                           group.by = 'orig.condition', 
                                                           subset.ident = clusterid,
                                                           test.use = DEG_TEST_USE)
  
  # Remove the possible infinite Log2FC values
  inf_log2fc_index = which( all_current_sample_cluster_markergenes_df$avg_log2FC == "Inf" | all_current_sample_cluster_markergenes_df$avg_log2FC == "-Inf")
  if( length( inf_log2fc_index) > 0){
    all_current_sample_cluster_markergenes_df = all_current_sample_cluster_markergenes_df[ -inf_log2fc_index,]
  }
  
  # Select genes with low enough p-value
  current_sample_cluster_markergenes_df = all_current_sample_cluster_markergenes_df[ which( all_current_sample_cluster_markergenes_df$p_val_adj <= MARKER_GENES_ALPHA_THRESHOLD), ]
  
  if( nrow( current_sample_cluster_markergenes_df) > 0){
    
    current_sample_cluster_markergenes_set = row.names( current_sample_cluster_markergenes_df)
    
    # Plot the volcano plot of the DEG
    print(
      EnhancedVolcano( toptable = all_current_sample_cluster_markergenes_df,
                       lab = row.names( all_current_sample_cluster_markergenes_df),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       selectLab = row.names( current_sample_cluster_markergenes_df),
                       xlab = bquote(~Log[2]~ 'fold change'),
                       pCutoff = 0.05,
                       FCcutoff = 1.2,
                       legendPosition = "right",
                       legendLabSize = 8,
                       drawConnectors = TRUE,
                       widthConnectors = 0.75,
                       title = paste( "DEG beetween", CONDITION_NAME_PTEN, "and", CONDITION_NAME_CTRL),
                       subtitle = paste( "Cluster", clusterid, ": PTEN=", nb_cell_PTEN_in_cluster, "cells / CTRL=", nb_cell_CTRL_in_cluster, "cells")
      )
    )
    
    current_sample_cluster_markergenes_df$gene = row.names( current_sample_cluster_markergenes_df)
    current_sample_cluster_markergenes_df$cluster = clusterid
    row.names( current_sample_cluster_markergenes_df) = NULL
    
    sample_cluster_markergenes_df = rbind( sample_cluster_markergenes_df, current_sample_cluster_markergenes_df)
    
  }else{
    current_sample_cluster_markergenes_set = NULL
    cat("<BR><b>No marker genes found for cluster", clusterid, " at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
  
  # Add a sub-tab with the enrichment of the positive markers in the GO BP terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### GO Biological Process enrichment in DEG")
  cat( "\n \n")
  
  # If some markers genes have been found, compute enrichment in GO BP
  if( ! is.null( current_sample_cluster_markergenes_set)){
    
    # Compute enrichments
    positive_ego <- enrichGO(gene          = current_sample_cluster_markergenes_set,
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
        
        # Save result to file
        positive_ego_df = data.frame( positive_ego)
        positive_ego_df$cluster = clusterid
        write.table( positive_ego_df, col.names = TRUE, row.names = FALSE, sep=";",
                     file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "GOBPEnrichment_Cluster", clusterid, ".csv"))
        )
        
        # Plot the results
        positive_dotplot = dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of marker genes in BP GO terms")
        ggsave( positive_dotplot, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "GOBPEnrichment_Cluster", clusterid, ".svg" )))
        print( positive_dotplot)
        
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
    cat("<BR><b>No marker genes found for cluster", clusterid, " at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
  
  
  # Add a sub-tab with the enrichment of the positive markers in the GO MF terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### GO Molecular Function enrichment in DEG")
  cat( "\n \n")
  
  # If some markers genes have been found, compute enrichment in GO MF
  if( ! is.null( current_sample_cluster_markergenes_set)){
    
    # Compute enrichments
    positive_ego <- enrichGO(gene          = current_sample_cluster_markergenes_set,
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
        
        positive_dotplot = dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of marker genes in MF GO terms")
        ggsave( positive_dotplot, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "GOMFEnrichment_Cluster", clusterid, ".svg" )))
        print( positive_dotplot)
        
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
  cat( "\n##### KEGG pathway enrichment in DEG")
  cat( "\n \n")
  
  # If some markers genes have been found, compute enrichment in GO MF
  if( ! is.null( current_sample_cluster_markergenes_set)){
    
    # Get the entrezid of the Seurat object genes
    universe_entrez_id = AnnotationDbi::mapIds( org.Mm.eg.db, rownames( merge_seurat_object@assays$RNA), keytype="SYMBOL", column="ENTREZID")
    
    # Convert the marker gene symbols into entrezid
    markers_entrez_id = AnnotationDbi::mapIds( org.Mm.eg.db, current_sample_cluster_markergenes_set, keytype="SYMBOL", column="ENTREZID")
    
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
        kegg_dotplot = dotplot( positive_kegg) + ggtitle("Enrichment of positive markers in chosen KEGG pathways")
        ggsave( kegg_dotplot, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "ChosenKEGGEnrichment_Cluster", clusterid, ".svg" )))
        print( kegg_dotplot)   
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
        kegg_dotplot = dotplot( positive_kegg) + ggtitle("Enrichment of positive markers in chosen KEGG pathways")
        ggsave( kegg_dotplot, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "ChosenKEGGEnrichment_Cluster", clusterid, ".svg" )))
        print( kegg_dotplot)   
      }else{
        cat("<BR>No enrichment found")
      }
    }else{
      cat("<BR>Enrichment analysis got an issue.")
      cat("<BR>markers_entrez_id = ", paste( markers_entrez_id, collapse = ";"))
    }
    
    cat(" \n \n"); # Required for '.tabset'
  }else{
    cat("<BR><b>No marker genes found for cluster", clusterid, " at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
  
}

all_sample_cluster_markergenes = unique( sample_cluster_markergenes_df$gene)

write.csv( sample_cluster_markergenes_df,
           file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "clusterDEG_", CONDITION_NAME_PTEN, "vs", CONDITION_NAME_CTRL, "_", DEG_TEST_USE, ".csv")),
           row.names = sample_cluster_markergenes_df$gene)



