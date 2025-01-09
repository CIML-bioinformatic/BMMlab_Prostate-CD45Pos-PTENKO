# ##############################################################
# This script aims to analyse the composition of clusters
# across condition and time
# ##############################################################

# .....................................................
## @knitr composition_analysis
# .....................................................

# See https://github.com/stemangiola/sccomp


# .....................................................
# BY CONDITION
# .....................................................


if( length( unique( merge_seurat_object$orig.condition)) > 1){
  cat("\n \n")
  cat("### By condition {.tabset .tabset-fade}")
  cat("\n \n")
  
  # Plot the cell in UMAP colored by clusters and split by condition
  print(
    DimPlot( merge_seurat_object, 
             group.by = "seurat_clusters", 
             split.by = "orig.condition",
             cols = CLUSTER_COLOR_PANEL[ 1:length( unique( merge_seurat_object$seurat_clusters))],
             label = TRUE)
  )

  cat("\n \n")
  cat("#### Barplot")
  cat("\n \n")
  
  # Look at the number of cells per condition and cluster
  Idents( merge_seurat_object) = "seurat_clusters"
  cells_clusterid = Idents( merge_seurat_object)
  Idents( merge_seurat_object) = "orig.condition"
  cells_condition = Idents( merge_seurat_object)
  
  # Look at the distribution of condition across clusters
  condition_vs_cluster_table = table( cells_clusterid, cells_condition)
  condition_vs_cluster_chisq_test = chisq.test( condition_vs_cluster_table)

  condition_vs_cluster_df = as.data.frame.table( condition_vs_cluster_table)
  print(
    ggplot( condition_vs_cluster_df, aes(fill=cells_clusterid, y=Freq, x=cells_condition)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual( values = CLUSTER_COLOR_PANEL[ 1 : length( levels( cells_clusterid))])
  )

  cat("\n \n")
  cat("#### Composition analysis")
  cat("\n \n")
  
  # Analyze the data composition by clusters against conditions
  sccomp_result_var = merge_seurat_object |>
                        sccomp_estimate( 
                          formula_composition = ~ orig.condition, 
                          formula_variability = ~ orig.condition,
                          .sample =  orig.ident, 
                          .cell_group = seurat_clusters,
                          bimodal_mean_variability_association = TRUE,
                          cores = 1 
                        ) |> 
                        sccomp_remove_outliers() |> 
                        sccomp_test(test_composition_above_logit_fold_change = 0.2)
  
  # Print the result table
  print(
    sccomp_result_var %>%kable() %>% 
      kable_styling( full_width = FALSE) %>%
      kable_paper("hover", full_width = T) %>%
      scroll_box( height = "500px")
  )
  
  
  # Get the plots from the results
  plots = plot( sccomp_result_var) 
  
  # Plot the posterior predictive check
  print( plots$boxplot)
  
  # Plot the 1D significance plots
  print( plots$credible_intervals_1D)
  
  # Save result to RDS and plot to SVG
  saveRDS( object = sccomp_result_var, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "CompositionAnalysis_ByCondition_ScCompObject.RDS")))
  
  ggsave( filename = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "CompositionAnalysis_ByCondition_IntervalPlot.svg")),
          plot = plots$credible_intervals_1D)
}


# .....................................................
# BY DATE
# .....................................................

if( length( unique( merge_seurat_object$orig.date)) > 1){
  cat("\n \n")
  cat("### By date {.tabset .tabset-fade}")
  cat("\n \n")
  
  # Plot the cell in UMAP colored by clusters and split by date
  print( 
    DimPlot( merge_seurat_object, 
             group.by = "seurat_clusters", 
             split.by = "orig.date",
             cols = CLUSTER_COLOR_PANEL[ 1:length( unique( merge_seurat_object$seurat_clusters))],
             label = TRUE)
  )
  
  cat("\n \n")
  cat("#### Barplot")
  cat("\n \n")
  
  # Look at the number of cells per date and cluster
  Idents( merge_seurat_object) = "seurat_clusters"
  cells_clusterid = Idents( merge_seurat_object)
  Idents( merge_seurat_object) = "orig.date"
  cells_date = Idents( merge_seurat_object)
  
  # Look at the distribution of date across clusters
  date_vs_cluster_table = table( cells_clusterid, cells_date)
  date_vs_cluster_chisq_test = chisq.test( date_vs_cluster_table)

  date_vs_cluster_df = as.data.frame.table( date_vs_cluster_table)
  print(
    ggplot( date_vs_cluster_df, aes(fill=cells_clusterid, y=Freq, x=cells_date)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual( values = CLUSTER_COLOR_PANEL[ 1 : length( levels( cells_clusterid))])
  )
  
  cat("\n \n")
  cat("#### Composition analysis")
  cat("\n \n")
  
  # Analyze the data composition by clusters against dates
  sccomp_result_var = merge_seurat_object |>
    sccomp_estimate( 
      formula_composition = ~ orig.date, 
      formula_variability = ~ orig.date,
      .sample =  orig.ident, 
      .cell_group = seurat_clusters,
      bimodal_mean_variability_association = TRUE,
      cores = 1 
    ) |> 
    sccomp_remove_outliers() |> 
    sccomp_test(test_composition_above_logit_fold_change = 0.2)
  
  # Print the result table
  print(
    sccomp_result_var %>%kable() %>% 
      kable_styling( full_width = FALSE) %>%
      kable_paper("hover", full_width = T) %>%
      scroll_box( height = "500px")
  )
  
  
  # Get the plots from the results
  plots = plot( sccomp_result_var) 
  
  # Plot the posterior predictive check
  print( plots$boxplot)
  
  # Plot the 1D significance plots
  print( plots$credible_intervals_1D)
  
  # Save result to RDS and plot to SVG
  saveRDS( object = sccomp_result_var, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "CompositionAnalysis_ByDate_ScCompObject.RDS")))
  
  ggsave( filename = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "CompositionAnalysis_ByDate_IntervalPlot.svg")),
          plot = plots$credible_intervals_1D)
}
