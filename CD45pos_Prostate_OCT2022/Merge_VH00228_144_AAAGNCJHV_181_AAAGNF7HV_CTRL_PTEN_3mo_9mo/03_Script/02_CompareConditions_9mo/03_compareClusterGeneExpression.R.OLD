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
  cat("#### Cluster", clusterid, "\n")
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
                       legendLabSize = 3,
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
    cat("<BR><b>No marker genes found for cluster", clusterid, " at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
}


all_sample_cluster_markergenes = unique( sample_cluster_markergenes_df$gene)

write.csv( sample_cluster_markergenes_df,
           file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "clusterDEG_", CONDITION_NAME_PTEN, "vs", CONDITION_NAME_CTRL, "_", DEG_TEST_USE, ".csv")),
           row.names = sample_cluster_markergenes_df$gene)



