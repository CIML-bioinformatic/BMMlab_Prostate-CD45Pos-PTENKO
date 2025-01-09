# ###################################################
# This script analyses the expression of gene modules
# across clusters and samples
# ###################################################


# Compare module scores between samples
# ###############################################

## @knitr compare_module_score

cat("\n \n")
cat("### Module score", CONDITION_NAME_PTEN, "vs", CONDITION_NAME_CTRL, "in clusters {.tabset .tabset-fade}\n")

PTEN_cells = names( which( merge_seurat_object$orig.condition == CONDITION_NAME_PTEN))
CTRL_cells = names( which( merge_seurat_object$orig.condition == CONDITION_NAME_CTRL))

for( current_module in names( MODULES_GENES)){
  
  # Check if the gene module has its score in the Seurat Object. If not, bypass it
  module_score_name = paste0( current_module, "1")
  if( ! module_score_name %in% names( merge_seurat_object@meta.data)){
    next
  }
  
  # Add the tab for this module
  cat("\n \n")
  cat("#### ",current_module)
  cat("\n \n")
  
  
  # Get the data of module score across clusters and conditions
  module_df = merge_seurat_object@meta.data[ , c( "seurat_clusters", "orig.condition", module_score_name)]
  module_df$orig.condition = factor( module_df$orig.condition, levels = c( CONDITION_NAME_CTRL, CONDITION_NAME_PTEN))
  
  ## Plot the module scoe across cluster and conditions
  print(
    ggplot( module_df, 
            aes_string( x="seurat_clusters", y=module_score_name, fill = "orig.condition", color = "orig.condition")) +
      geom_violin( position = "dodge") + 
      stat_summary( fun = mean, geom = "point", color = "black", position = position_dodge( width = 0.9)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", size = 0.5, width = 0.3, position = position_dodge( width = 0.9)) +
      scale_fill_manual( values = CONDITION_COLOR[ c( CONDITION_NAME_CTRL, CONDITION_NAME_PTEN)]) +
      scale_color_manual( values = CONDITION_COLOR[ c( CONDITION_NAME_CTRL, CONDITION_NAME_PTEN)]) +
      geom_hline( yintercept=0, linetype="dashed", color = "black", linewidth = 1) + 
      theme_minimal() + ggtitle( paste( "Module scores of", current_module, "module\nby cluster and condition with mean and standard deviation"))
  )
  
  # Compute a global Kruskall-Wallis statistics on the global data
  module_df$orig.condition_x_seurat_clusters = paste( module_df$orig.condition, module_df$seurat_clusters, sep = "-")
  k_test = kruskal.test( formula = as.formula( paste( module_score_name, " ~ orig.condition_x_seurat_clusters")), 
                         data = module_df )
  
  cat("<BR>Kruskal-Wallis chi-squared =", k_test$statistic, ", df =", k_test$parameter[[ "df"]] , ", p-value = ", k_test$p.value, "<BR>")
  
  # Execute the post-hoc analysis using wilcoxon test comparing for each cluster the two conditions
  module_stats_df = data.frame()
  for( clusterid in levels( merge_seurat_object@meta.data$seurat_clusters)){
    data_cluster = merge_seurat_object@meta.data[ which( merge_seurat_object@meta.data$seurat_clusters == clusterid), c( "seurat_clusters", "orig.condition", module_score_name)]
    if( sum( data_cluster[ , module_score_name]) == 0){
      cat("<BR>Note : No statistic for cluster", clusterid, "because all values are zero")
    }else{
      nb_CTRL_cells = length( which( data_cluster$orig.condition == CONDITION_NAME_CTRL))
      nb_PTEN_cells = length( which( data_cluster$orig.condition == CONDITION_NAME_PTEN))
      if( nb_CTRL_cells >= 10 && nb_PTEN_cells >= 10){
        stats_df = wilcox_test( data = data_cluster, formula = as.formula( paste( module_score_name, "~ orig.condition")))
        stats_df$seurat_clusters = clusterid
        module_stats_df = rbind( module_stats_df, stats_df)        
      }else{
        cat("<BR>To few cells in cluster", clusterid, "to perform analysis")
      }
    }
    
  }
  
  if( nrow( module_stats_df) > 0){
    # Adjust the p-value for multi-testing with BH method to get FDR
    module_stats_df$FDR = signif( p.adjust( module_stats_df$p, method = "BH"), 3)
    module_stats_df$signif = sapply( module_stats_df$FDR, function( fdr){ if( fdr < 0.01) "**" else if( fdr < 0.05) "*" else "" })
    module_stats_df = module_stats_df[ , c( "seurat_clusters", "group1", "group2", "n1", "n2", "p", "FDR", "signif")]
    
    # Display the statistics result table
    print(
      module_stats_df %>% 
        kable( caption = paste( current_module, " : Comparison between conditions by clusters with Wilcoxon test")) %>% 
        kable_styling( full_width = FALSE)
    )
  }else{
    cat("<BR>No statistic summary table to display since no statistics were computed.")
  }
  
}
