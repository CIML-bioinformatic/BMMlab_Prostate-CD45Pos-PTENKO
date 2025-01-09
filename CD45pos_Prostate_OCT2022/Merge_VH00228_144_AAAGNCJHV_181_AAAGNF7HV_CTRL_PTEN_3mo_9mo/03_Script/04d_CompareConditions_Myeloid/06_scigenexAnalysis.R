# #######################################################################
# This script aims to apply SciGeneX analysis to evaluate the
# interest of this approach on the datasets
# #######################################################################

# SciGeneX analysis
# ###############################################

## @knitr scigenex_analysis

# Set the verbosity level to 0 in order to avoid messages in the HTM Lreport.
scigenex::set_verbosity(0)

# Select informative genes
scigenex_analysis = select_genes(  merge_seurat_object,
                                   distance = "pearson",
                                   row_sum=5)

# Construct and partition the graph
scigenex_analysis <- gene_clustering( scigenex_analysis,
                                      inflation = 2,
                                      threads = 4)

# Display the heatmap of gene clusters
cat("\n### Gene module before filtering\n")
scigenex_analysis <- top_genes( scigenex_analysis)

Idents( merge_seurat_object) = "seurat_clusters"
plot_heatmap( scigenex_analysis,
              interactive = FALSE,
              cell_clusters = Idents( merge_seurat_object))

# Display information about the gene modules
plot_cluster_stats( cluster_stats(scigenex_analysis))

# Filter gene modules
##... based on number of genes within modules
scigenex_analysis <- filter_cluster_size(scigenex_analysis, min_cluster_size = 5)

#... based on standard deviation
scigenex_analysis <- filter_cluster_sd(scigenex_analysis, min_sd = 0.4)

cat("\n### Gene module after filtering\n")

scigenex_analysis <- top_genes( scigenex_analysis)

Idents( merge_seurat_object) = "seurat_clusters"
plot_heatmap( scigenex_analysis,
              interactive = FALSE,
              cell_clusters = Idents( merge_seurat_object))

cluster_condition = paste( unlist( merge_seurat_object[[ "seurat_clusters"]], use.names = TRUE),
                           unlist( merge_seurat_object[[ "orig.condition"]], use.names = TRUE), sep="_")
names( cluster_condition) = Cells( merge_seurat_object)
cluster_condition = factor( cluster_condition,
                            levels = sort( unique( cluster_condition)))

merge_seurat_object = AddMetaData( merge_seurat_object, 
                                   metadata = cluster_condition,
                                   col.name = "Cluster_Condition")


Idents( merge_seurat_object) = "Cluster_Condition"
plot_heatmap( scigenex_analysis,
              interactive = FALSE,
              cell_clusters = Idents( merge_seurat_object))

# Display information about the gene modules
plot_cluster_stats( cluster_stats(scigenex_analysis))

# Get the list of top genes and save them to file
module_max_size = max( sapply( scigenex_analysis@top_genes, length))
gene_list_set = c()
for( current_module in names( scigenex_analysis@top_genes)){
  gene_list = scigenex_analysis@top_genes[[ current_module]]
  if( length( gene_list) < module_max_size){
    gene_list = c( gene_list, rep( "", module_max_size - length( gene_list)))
  }
  gene_list_set = c( gene_list_set, paste( gene_list, collapse=";"))
}

fileConn<-file( file.path( PATH_ANALYSIS_EXTRA_OUTPUT, "SciGeneX_TopGenes"))
writeLines( gene_list_set, fileConn)
close(fileConn)
