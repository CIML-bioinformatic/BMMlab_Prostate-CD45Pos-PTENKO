# ###################################################
# This script analyses the expression of gene modules
# across clusters and samples
# ###################################################


# Compare module scores between samples
# ###############################################

## @knitr compare_module_score

cat("\n \n")
cat("### Module score", DATE_NAME_9MO, "vs", DATE_NAME_3MO, "in clusters {.tabset .tabset-fade}\n")

date9mo_cells = names( which( merge_seurat_object$orig.date == DATE_NAME_9MO))
date3mo_cells = names( which( merge_seurat_object$orig.date == DATE_NAME_3MO))

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
  
  
  # Get the data of module score across clusters and dates
  module_df = merge_seurat_object@meta.data[ , c( "seurat_clusters", "orig.date", module_score_name)]
  module_df$orig.date = factor( module_df$orig.date, levels = c( DATE_NAME_3MO, DATE_NAME_9MO))
  
  ## Plot the module scoe across cluster and dates
  print(
    ggplot( module_df, 
            aes_string( x="seurat_clusters", y=module_score_name, fill = "orig.date", color = "orig.date")) +
      geom_violin( position = "dodge") + 
      stat_summary( fun = mean, geom = "point", color = "black", position = position_dodge( width = 0.9)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", size = 0.5, width = 0.3, position = position_dodge( width = 0.9)) +
      scale_fill_manual( values = DATE_COLOR[ c( DATE_NAME_3MO, DATE_NAME_9MO)]) +
      scale_color_manual( values = DATE_COLOR[ c( DATE_NAME_3MO, DATE_NAME_9MO)]) +
      geom_hline( yintercept=0, linetype="dashed", color = "black", linewidth = 1) + 
      theme_minimal() + ggtitle( paste( "Module scores of", current_module, "module\nby cluster and date with mean and standard deviation"))
  )
  
  # Compute a global Kruskall-Wallis statistics on the global data
  module_df$orig.date_x_seurat_clusters = paste( module_df$orig.date, module_df$seurat_clusters, sep = "-")
  k_test = kruskal.test( formula = as.formula( paste( module_score_name, " ~ orig.date_x_seurat_clusters")), 
                         data = module_df )
  
  cat("<BR><b>Kruskal-Wallis chi-squared =", k_test$statistic, ", df =", k_test$parameter[[ "df"]] , ", p-value = ", k_test$p.value, "</b><BR>")
  
  # Execute the post-hoc analysis using wilcoxon test comparing for each cluster the two dates
  module_stats_df = data.frame()
  for( clusterid in levels( merge_seurat_object@meta.data$seurat_clusters)){
    data_cluster = merge_seurat_object@meta.data[ which( merge_seurat_object@meta.data$seurat_clusters == clusterid), c( "seurat_clusters", "orig.date", module_score_name)]
    cluster_3mo_cells_nb = length( which( data_cluster$seurat_clusters == clusterid & data_cluster$orig.date == DATE_NAME_3MO))
    cluster_9mo_cells_nb = length( which( data_cluster$seurat_clusters == clusterid & data_cluster$orig.date == DATE_NAME_9MO))
    if( cluster_3mo_cells_nb > 10 && cluster_9mo_cells_nb > 10){
      stats_df = wilcox_test( data = data_cluster, formula = as.formula( paste( module_score_name, "~ orig.date")))
      stats_df$seurat_clusters = clusterid
      module_stats_df = rbind( module_stats_df,
                              stats_df)
    }else{
      cat("<BR><b>Cluster", clusterid, ": no statistics because de nb of cells is too low (", cluster_3mo_cells_nb,"for 3mo cells,", cluster_9mo_cells_nb, "for 9mo cells)</b>")
    }
                            
  }
  
  # Adjust the p-value for multi-testing with BH method to get FDR
  module_stats_df$FDR = signif( p.adjust( module_stats_df$p, method = "BH"), 3)
  module_stats_df$signif = sapply( module_stats_df$FDR, function( fdr){ if( fdr < 0.01) "**" else if( fdr < 0.05) "*" else "" })
  module_stats_df = module_stats_df[ , c( "seurat_clusters", "group1", "group2", "n1", "n2", "p", "FDR", "signif")]
  
  # Display the statistics result table
  print(
    module_stats_df %>% 
      kable( caption = paste( current_module, " : Comparison between dates by clusters with Wilcoxon test")) %>% 
      kable_styling( full_width = FALSE)
  )
  
}
