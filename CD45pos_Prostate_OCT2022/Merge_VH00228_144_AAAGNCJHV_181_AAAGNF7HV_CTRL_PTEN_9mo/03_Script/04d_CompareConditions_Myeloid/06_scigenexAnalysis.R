# #######################################################################
# This script aims to apply SciGeneX analysis to evaluate the
# interest of this approach on the datasets
# #######################################################################

# SciGeneX analysis
# ###############################################

## @knitr scigenex_analysis

# .................................................................................
# Launch the SciGeneX analysis
# .................................................................................

cat("\n \n")
cat("### Gene module analysis", module_id)
cat("\n \n")

# Set the verbosity level to 0 in order to avoid messages in the HTM Lreport.
scigenex::set_verbosity(0)

# Select informative genes
scigenex_analysis = select_genes(  merge_seurat_object,
                                   distance = "pearson",
                                   row_sum=5)

## Construct and partition the graph
scigenex_analysis_1_5 <- gene_clustering( scigenex_analysis,
                                      inflation = 1.5,
                                      threads = 4)

# .................................................................................
# Look at the gene modules in heatmap before and after filtering
# .................................................................................

cat("\n \n")
cat("#### Heatmap before filtering")
cat("\n \n")

# Display the heatmap of gene clusters
scigenex_analysis <- top_genes( scigenex_analysis)
Idents( merge_seurat_object) = "seurat_clusters"
plot_heatmap( scigenex_analysis,
              interactive = FALSE,
              cell_clusters = Idents( merge_seurat_object))

# Display information about the gene modules
plot_cluster_stats( cluster_stats(scigenex_analysis))

# Generate html report (contains interactive heatmap)
# unfiltered_report_file_path = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, "SciGeneX_interative_heatmap_notFiltered.html")
# cluster_set_report(clusterset_object = scigenex_analysis,
#                    seurat_object = merge_seurat_object,
#                    file_path = unfiltered_report_file_path,
#                    force = TRUE,
#                    verbosity = FALSE)
# cat("<BR><b>Report on unfiltered gene modules has been saved to:\n</b>", unfiltered_report_file_path, "<BR>")

cat("\n \n")
cat("#### Heatmap after filtering")
cat("\n \n")

# Filter gene modules
##... based on number of genes within modules
scigenex_analysis <- filter_cluster_size(scigenex_analysis, min_cluster_size = 20)

#... based on standard deviation
scigenex_analysis <- filter_cluster_sd(scigenex_analysis, min_sd = 0.4)

plot_heatmap( scigenex_analysis,
              interactive = FALSE,
              cell_clusters = Idents( merge_seurat_object))

# .................................................................................
# Get the top genes of each module and save them to file
# .................................................................................

# Get the list of all genes and save them to file
gene_module_list = c()
for( current_module in names( scigenex_analysis@gene_clusters)){
  gene_list = get_genes(scigenex_analysis, cluster = current_module)
  gene_module_list = c( gene_module_list, paste0( "Module_", current_module, ";",paste0(  gene_list, collapse=";")))
}

all_gene_file_path = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, "SciGeneX_AllGenes.csv")
fileConn<-file( all_gene_file_path)
writeLines( gene_module_list, fileConn)
close(fileConn)

cat("<BR>The list of top genes per module has been saved to file:<BR>", top_gene_file_path, "<BR>")

# Get the list of top genes and save them to file
gene_module_list = c()
for( current_module in names( scigenex_analysis@top_genes)){
  gene_list = scigenex_analysis@top_genes[[ current_module]]
  gene_module_list = c( gene_module_list, paste0( "Module_", current_module, ";",paste0(  gene_list, collapse=";")))
}

top_gene_file_path = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, "SciGeneX_TopGenes.csv")
fileConn<-file( top_gene_file_path)
writeLines( gene_module_list, fileConn)
close(fileConn)

cat("<BR>The list of top genes per module has been saved to file:<BR>", top_gene_file_path, "<BR>")

# .................................................................................
# Look at the module score of each gene module
# .................................................................................
 
cat("\n \n")
cat("### Gene module scores (all genes)")
cat("\n \n")

for( module_id in names( scigenex_analysis@gene_clusters)){
  
  
  cat("\n \n")
  cat("#### Module", module_id)
  cat("\n \n")
  module_name = paste0( "SciModule_", module_id, "_")
  
  module_gene_set = get_genes(scigenex_analysis, cluster = module_id)
  merge_seurat_object = AddModuleScore( object = merge_seurat_object, 
                                        features = list( module_gene_set),
                                        name = module_name)
  print( FeaturePlot( merge_seurat_object, features = paste0( module_name, "1")) +
           ggtitle( paste( "Score of Module", module_id, "by cluster"))
  )
  
  print( VlnPlot( merge_seurat_object, features = paste0( module_name, "1")) + 
    geom_hline(yintercept = 0, linetype = 2, color = "red", linewidth = 1) +
    ggtitle( paste( "Score of Module", module_id, "by cluster"))
  )
  cat("\n \n")
}

cat("\n \n")

cat("\n \n")
cat("### Gene module scores (top genes)")
cat("\n \n")

for( module_id in names( scigenex_analysis@top_genes)){
  
  
  cat("\n \n")
  cat("#### Module", module_id)
  cat("\n \n")
  module_name = paste0( "SciModule_", module_id, "_TopGenes_")
  
  module_gene_set = scigenex_analysis@top_genes[[ module_id]]
  merge_seurat_object = AddModuleScore( object = merge_seurat_object, 
                                        features = list( module_gene_set),
                                        name = module_name)
  print( FeaturePlot( merge_seurat_object, features = paste0( module_name, "1")) +
           ggtitle( paste( "Score of Module", module_id, "(top genes) by cluster"))
         )
  
  print( VlnPlot( merge_seurat_object, features = paste0( module_name, "1")) + 
    geom_hline(yintercept = 0, linetype = 2, color = "red", linewidth = 1) +
    ggtitle( paste( "Score of Module", module_id, "(top genes) by cluster"))
  )
  cat("\n \n")
}

cat("\n \n")





