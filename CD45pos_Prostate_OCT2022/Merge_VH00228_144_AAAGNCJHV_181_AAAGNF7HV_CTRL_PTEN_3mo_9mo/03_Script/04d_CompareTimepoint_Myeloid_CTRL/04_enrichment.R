# ##################################################
# This script reads data from previous analysis
# ##################################################




# Compare expression of gene in cluster between samples
# -------------------------------------------------------

## @knitr compare_cluster_expression

cat("\n \n")
cat("### DEG between KO and WT in clusters {.tabset .tabset-fade}\n")

# Detect the DEG between samples in each cluster
# ...............................................
sample_cluster_markergenes_df = data.frame()

Idents( merge_seurat_object) = "seurat_clusters"
for( clusterid in levels( Idents( merge_seurat_object))){
  
  cat("\n \n")
  cat("#### Cluster", clusterid, "\n")
  
  all_current_sample_cluster_markergenes_df = FindMarkers( merge_seurat_object, 
                                                               ident.1 = SAMPLE_KO, 
                                                               group.by = 'HTO_classification', 
                                                               subset.ident = clusterid,
                                                               test.use = DEG_TEST_USE)
  
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
                       drawConnectors = TRUE,
                       widthConnectors = 0.75,
                       title = paste( "DEG beetween", SAMPLE_KO, "and", SAMPLE_WT),
                       subtitle = paste( "Cluster", clusterid)
      )
    )
    
    
    current_sample_cluster_markergenes_df$gene = row.names( current_sample_cluster_markergenes_df)
    current_sample_cluster_markergenes_df$cluster = clusterid
    row.names( current_sample_cluster_markergenes_df) = NULL
    
    sample_cluster_markergenes_df = rbind( sample_cluster_markergenes_df, current_sample_cluster_markergenes_df)

    print(
      current_sample_cluster_markergenes_df %>% kable( caption = "List of DEG genes by cluster") %>% kable_styling()
    )
    
  }else{
    cat("<BR><b>No marker genes found for cluster", clusterid, " at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
}


all_sample_cluster_markergenes = unique( sample_cluster_markergenes_df$gene)

write.csv( sample_cluster_markergenes_df,
           file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "clusterDEG_KOvsWT_", DEG_TEST_USE, ".csv")),
           row.names = sample_cluster_markergenes_df$gene)


