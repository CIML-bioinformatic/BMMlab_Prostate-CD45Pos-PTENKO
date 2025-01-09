# ####################################################
# This script aim to use Monocle3 pipeline
# to analyse the cell trajectory
# ####################################################

## @knitr estimate_trajectory


# Trajectory analysis

# A subtab for the clustering
# ...............................................................................
cat( "\n \n")
cat( "\n### Clusters")
cat( "\n \n")

# -- Cluster the cells
cds_sc10X_seurat = cluster_cells( cds_sc10X_seurat, resolution = CLUSTER_RESOLUTION)

# -- Plot the clusters on UMAP
p1 = plot_cells( cds_sc10X_seurat, 
            color_cells_by = "cluster", 
            group_label_size = 10,
            show_trajectory_graph = FALSE, 
            cell_size = 1) + ggtitle( paste( "Cluster distribution\n for resolution", CLUSTER_RESOLUTION))

print( p1)
cat( "<BR>")


# A subtab for the trajectory
# ...............................................................................
cat( "\n \n")
cat( "\n### Trajectory Graph")
cat( "\n \n")

# -- Learn the trajectory from clusters and expression data
cds_sc10X_seurat = learn_graph( cds_sc10X_seurat, use_partition = FALSE, verbose = FALSE)

# -- Plot the trajectory on UMAP
trajectory_plot = plot_cells( cds_sc10X_seurat,
                      color_cells_by = "cluster",
                      label_groups_by_cluster=FALSE,
                      label_leaves=FALSE,
                      label_branch_points=FALSE, 
                      label_principal_points = TRUE,
                      graph_label_size = 4,
                      cell_size = 0.5)  +
                    ggtitle( paste( "Trajectory for\n cluster resolution", CLUSTER_RESOLUTION))
print( trajectory_plot)

ggsave( filename = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "MonocleTrajectories.svg")),
        trajectory_plot,
        bg = "white")

cat( "<BR>")


# ...............................................................................
# ...............................................................................
# For some resolution values, get some extra details
# ...............................................................................
# ...............................................................................

cat( "\n \n")
cat( "\n### Trajectories detail  {.tabset}")
cat( "\n \n")

# Identify the nereast cell of each principal vertex
closest_vertex = cds_sc10X_seurat@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
closest_to = list()
for( vertex in sort( unique( closest_vertex[ ,1]))){
  closest_to[[ paste0( "Y_", vertex)]] = row.names(closest_vertex)[ which( closest_vertex == vertex)]
}  

# Identify the vertex of degree 1 in the principal graph (source or destination)
vertex_degrees <- degree( cds_sc10X_seurat@principal_graph$UMAP)
vertex_degree_1 = names( vertex_degrees[ vertex_degrees == 1])
pairs_vertex_degree_1 = combn( vertex_degree_1, 2)

# Look at the genes to monitor
#filtered_monitored_genes = intersect( MONITORED_GENES, row.names( sc10x.seurat@assays$RNA@data))

# Look at the cells for each condition
Idents( sc10x.seurat) = "orig.date"
date_3mo_cells = names( Idents( sc10x.seurat))[ which( Idents( sc10x.seurat) == DATE_NAME_3MO)]
date_9mo_cells = names( Idents( sc10x.seurat))[ which( Idents( sc10x.seurat) == DATE_NAME_9MO)]

# For each pair of source/destination node, analyse the expression of gene along the
# shortest path between the nodes
for( pair_index in 1:ncol( pairs_vertex_degree_1)){
  
  origin = pairs_vertex_degree_1[ 1, pair_index]
  destination = pairs_vertex_degree_1[ 2, pair_index]
  
  cat( "\n \n")
  cat( "\n#### Trajectory:", origin, "to", destination, " {.tabset}")
  cat( "\n \n")
  
  # Create a folder to store the plots on the details of this trajectory
  path_trajectory_details = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( "Trajectory_", origin, "_", destination))
  dir.create( path_trajectory_details, recursive = TRUE, showWarnings = FALSE)
  
  # Identify the shortest path between two principal vertex
  principal_graph = cds_sc10X_seurat@principal_graph@listData$UMAP
  shortest_path = shortest_paths( principal_graph, from = origin, to = destination, mode="all", output = "both", predecessors = TRUE)
  
  # For each principal vertex in path, retrieve the nearest cells
  path_cells = data.frame()
  for( vertex_index in 1:length( shortest_path$vpath[[1]])){
    vertex_name =  as.character( names( shortest_path$vpath[[1]][ vertex_index]))
    nearest_cells = closest_to[[vertex_name]]
    if( ! is.null(nearest_cells)){
      path_cells = rbind( path_cells, data.frame( principal.vertex.id = vertex_name, cell.id = nearest_cells))
    }
  }
  path_cells$principal.vertex.id = factor( path_cells$principal.vertex.id, levels = names( shortest_path$vpath[[1]]))
  
  # Compute the pseudotime of cells from the nearest cells of the origin principal vertex
  cds_sc10X_seurat <- order_cells( cds_sc10X_seurat, root_pr_nodes = origin)
  path_cells$pseudotime = cds_sc10X_seurat@principal_graph_aux@listData$UMAP$pseudotime[ path_cells$cell.id]

  # Provide the score of the modules of genes    
  module_score = sc10x.seurat@meta.data[ path_cells$cell.id, paste0( names( MODULES_GENES), "1")]
  module_score$cell.id = row.names( module_score)
  path_module_score = merge( path_cells, module_score, by = "cell.id")
  path_module_score = path_module_score[ order( path_module_score$pseudotime, decreasing = FALSE), ]
      
  # # Provide the score of the group of genes    
  # group_score = sc10x.seurat@meta.data[ path_cells$cell.id, paste0( names( GROUP_GENES), "1")]
  # group_score$cell.id = row.names( group_score)
  # path_group_score = merge( path_cells, group_score, by = "cell.id")
  # path_group_score = path_group_score[ order( path_group_score$pseudotime, decreasing = FALSE), ]
  
  # Inject the nearest principal vertex information as metadata of the Seurat object
  row.names( path_cells) = path_cells$cell.id
  meta_principal_vertex = rep(NA, length( Cells( sc10x.seurat)))
  names( meta_principal_vertex) = Cells( sc10x.seurat)
  meta_principal_vertex[ row.names( path_cells)] = as.character( path_cells$principal.vertex.id)
  meta_principal_vertex = factor( meta_principal_vertex, levels=names( shortest_path$vpath[[1]]))
  
  sc10x.seurat = AddMetaData( sc10x.seurat, metadata = meta_principal_vertex, col.name = "principal.vertex.id")
  Idents( sc10x.seurat) = "principal.vertex.id"
  print( DimPlot( sc10x.seurat) + 
           # theme( legend.position = "None")) + 
    ggtitle( "Nearest cells to principal vertices") + 
    scale_colour_viridis_d( na.value = "grey80")
  )
  
  # Inject the pseudotime information as metadata of the Seurat object
  meta_pseudotime = rep( NA, length( Cells( sc10x.seurat)))
  names( meta_pseudotime) = Cells( sc10x.seurat)
  meta_pseudotime[ row.names( path_cells)] = as.numeric( path_cells$pseudotime)
  
  sc10x.seurat = AddMetaData( sc10x.seurat, metadata = meta_pseudotime, col.name = "pseudotime")
  pseudotime_feature_plot =  FeaturePlot( sc10x.seurat, features = c( pseudotime = "pseudotime")) +
                                         scale_colour_viridis_c( option = "plasma", na.value = "grey80") +
                                         ggtitle( paste( "Cells pseudotime", origin, "to", destination))
  print( pseudotime_feature_plot)
  ggsave( filename = file.path( path_trajectory_details, "pseudotime_feature_plot.svg"),
          pseudotime_feature_plot,
          bg = "white")
  
  cat("<BR>")
  
  cat( "\n \n")
  cat( "\n##### Module gene expression  {.tabset}")
  cat( "\n \n")
  
  # For each monitored gene, look at its expression along the trajectory
  for( gene_module in names( MODULES_GENES)){
    
    cat( "\n \n")
    cat( "\n######", gene_module)
    cat( "\n \n")
    
    # Plot the computed statistics and expression of gene along path
    gene_module_pseudotime_plot = ggplot( data = path_module_score, aes_string(  x="pseudotime", y = paste0( gene_module, "1"))) + 
                                   geom_point( aes( color = pseudotime)) + geom_smooth( formula = 'y ~ x', method = "loess") +
                                   scale_colour_viridis_c( option = "plasma") +
                                   theme_minimal() + ggtitle( paste( "Expression of module", gene_module, "\nalong pseudotime for all cells"))
    
    print( gene_module_pseudotime_plot)
    ggsave( filename = file.path( path_trajectory_details, paste0( "modulescore_vs_pseudotime_", gene_module, ".svg")),
            gene_module_pseudotime_plot,
            bg = "white")
# 
#       print( ggplot(  data = NULL, aes_string(  x="pseudotime", y = paste0( gene_module, "1"))) + 
#                geom_point( data = path_module_score[ which( path_module_score$cell.id %in% date_3mo_cells), ], color = DATE_COLOR[[ SELECT]]) +
#                geom_smooth( data = path_module_score[ which( path_module_score$cell.id %in% date_3mo_cells), ], formula = 'y ~ x', method = "loess", color = DATE_COLOR[[ DATE_NAME_3MO]]) + 
#                geom_point( data = path_module_score[ which( path_module_score$cell.id %in% date_9mo_cells), ], color = DATE_COLOR[[ DATE_NAME_9MO]]) +
#                geom_smooth( data = path_module_score[ which( path_module_score$cell.id %in% date_9mo_cells), ], formula = 'y ~ x', method = "loess", color = DATE_COLOR[[ DATE_NAME_9MO]]) + 
#                theme_minimal() + ggtitle( paste( "Expression of module", gene_module, "\nalong pseudotime for 3mo and 9mo cells"))
#       )
    
    cat("<BR>")
 
  }
  
  # cat( "\n \n")
  # cat( "\n##### Chosen groups of genes score {.tabset}")
  # cat( "\n \n")
  # 
  # # For each group of genes, look at its score along the trajectory
  # for( groupid in names( GROUP_GENES)){
  #   
  #   cat( "\n \n")
  #   cat( "\n######", groupid)
  #   cat( "\n \n")
  #   
  #   # Plot the scores of group of gene along path
  #   print( ggplot( data = path_group_score, aes_string(  x="pseudotime", y = paste0( groupid, "1"))) + 
  #            geom_point( aes( color = pseudotime)) + geom_smooth( formula = 'y ~ x', method = "loess") +
  #            scale_colour_viridis_c( option = "plasma") +
  #            theme_minimal() + ggtitle( paste( "Score of", groupid, "along pseudotime for all cells"))
  #   )
  #   
  #   print( ggplot(  data = NULL, aes_string(  x="pseudotime", y = paste0( groupid, "1"))) + 
  #            geom_point( data = path_group_score[ which( path_group_score$cell.id %in% WT_cells), ], color = SAMPLE_COLOR[[ SAMPLE_WT]]) +
  #            geom_smooth( data = path_group_score[ which( path_group_score$cell.id %in% WT_cells), ], formula = 'y ~ x', method = "loess", color = SAMPLE_COLOR[[ SAMPLE_WT]]) + 
  #            geom_point( data = path_group_score[ which( path_group_score$cell.id %in% KO_cells), ], color = SAMPLE_COLOR[[ SAMPLE_KO]]) +
  #            geom_smooth( data = path_group_score[ which( path_group_score$cell.id %in% KO_cells), ], formula = 'y ~ x', method = "loess", color = SAMPLE_COLOR[[ SAMPLE_KO]]) + 
  #            theme_minimal() + ggtitle( paste( "Score of", groupid, "along pseudotime for WT and KO cells"))
  #   )
  #   
  #   cat("<BR>")
  #   
  # }
}

