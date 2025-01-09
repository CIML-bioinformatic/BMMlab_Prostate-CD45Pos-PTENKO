# #########################################################################
# This scripts aims to display the score of modules of genes in UMAP
# #########################################################################


# ##################################################

## @knitr heterogeneity_modules

# ##################################################

module_information_df = data.frame()

for( module_name in names( MODULES_GENES)){
  
  human_module = FALSE
  initial_length = length( MODULES_GENES[[ module_name]])
  
  # detect if genes are in upper case of not
  if( sum( toupper( MODULES_GENES[[ module_name]]) == MODULES_GENES[[ module_name]]) > (length( MODULES_GENES[[ module_name]]) / 2)){
    human_module = TRUE
    convertion = convertHumanGeneList( MODULES_GENES[[ module_name]])
    if( !is.null( convertion) && nrow( convertion) > 0){
      MODULES_GENES[[ module_name]] = convertion$symbol
    }else{
      MODULES_GENES[[ module_name]] = NULL
    }
  }
  
  if( module_name %in% names( MODULES_GENES)){
    final_length = length( MODULES_GENES[[ module_name]])
    module_information_df = rbind( module_information_df, data.frame( human.module = human_module,
                                                                      initial.length = initial_length,
                                                                      final.length = final_length,
                                                                      gene.list = paste( MODULES_GENES[[ module_name]], collapse=';'))
                                   )
  }
}
row.names( module_information_df) = names( MODULES_GENES)

print( module_information_df %>% kable( caption = "Information on gene modules") %>% kable_styling())

### Remove eventual NULL (empty) list elements from list of genes in modules
modulesGroupEmpty = sapply( MODULES_GENES, is.null);
if(any( modulesGroupEmpty)) warning( paste("Following module(s) of genes will be ignored because empty:", paste( names(modulesGroupEmpty)[modulesGroupEmpty], collapse=" - ")));
MODULES_GENES = MODULES_GENES[! modulesGroupEmpty];

# Check whether genes in MODULES_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
matchModulesGenes = match( ( unlist( MODULES_GENES)), ( rownames( GetAssayData( merge_seurat_object))));
modulesGenesNotFound = unique( unlist( MODULES_GENES)[is.na( matchModulesGenes)]);
if(any( is.na( matchModulesGenes))) warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:", paste( paste0("'", modulesGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MODULES_GENES = relist( rownames( GetAssayData( merge_seurat_object))[ matchModulesGenes ], skeleton = MODULES_GENES); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
MODULES_GENES = lapply( MODULES_GENES, na.omit);

# Just remind the warning for genes names not in object, or modules that were transfered to individual monitoring of genes
if(any( is.na( matchModulesGenes))){
  warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:",
                  paste( modulesGenesNotFound, collapse=" - ")));
}


# ##################################################

## @knitr heterogeneity_modules_scoring

# ##################################################

# Compute the score of the cells according to group of monitored genes
for( moduleName in names( MODULES_GENES)){
  message(moduleName);
  if( length( MODULES_GENES[[moduleName]]) == 0)
  {
    warning( paste0( "List of genes in module '", moduleName, "' is empty, ignoring..."));
    MODULES_GENES[[moduleName]] = NULL
  } else
  {
    merge_seurat_object <- AddModuleScore( object = merge_seurat_object,
                             features = MODULES_GENES[ moduleName], # Must be a list
                             ctrl = MODULES_CONTROL_SIZE,           #length(MODULES_GENES[[ moduleName]]),
                             name = moduleName,
                             seed = SEED);
  }
}



# ##################################################

## @knitr heterogeneity_modules_heatmap

# ##################################################

# DoHeatmap replaced by use of iHeatmapr after testing several options

# Get the matrix of module scores (not expression) for each cell and associated clusters from Seurat object
modulesScoreMat = t( as.matrix( merge_seurat_object[[ paste0(names( MODULES_GENES),1)]]));
clusterID = Idents( merge_seurat_object);

# Remove the extra numeric character added to modules names by Seurat
rownames( modulesScoreMat) = substr( rownames( modulesScoreMat), 1, nchar( rownames( modulesScoreMat))-1);

# Select genes in modules and reorder cells to group clusters together
clusterOrdering = order( clusterID);

modulesScoreMat = modulesScoreMat[, clusterOrdering];
clusterID = clusterID[ clusterOrdering];

# Prepare rows and columns annotation bars (module and cluster respectively)
#rowsAnnot = data.frame( Module = names( MODULES_GENES));
colsAnnot = data.frame( Cluster = clusterID);


# Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
originalRowNames = rownames( modulesScoreMat);
originalColNames = colnames( modulesScoreMat);
rownames( modulesScoreMat) = make.unique( originalRowNames);
colnames( modulesScoreMat) = make.unique( originalColNames);
#rownames( rowsAnnot) = rownames( modulesScoreMat);
rownames( colsAnnot) = colnames( modulesScoreMat);

# Plot the 'non-interactive' heatmap
pheatmap( modulesScoreMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          #          annotation_row = rowsAnnot,
          annotation_col = colsAnnot,
          labels_row = originalRowNames,
          annotation_colors = list( Cluster = clustersColor),
          show_colnames = FALSE);



# ##################################################

## @knitr heterogeneity_modules_expression_projection

# ##################################################

# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot scores of modules on dimensionality reduction figures
invisible( lapply( names(MODULES_GENES), function(featureName)
{
  print( FeaturePlot( merge_seurat_object, features = paste0(featureName, "1"), reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
           ggtitle( label = featureName) +
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none"));
}));

cat(" \n \n"); # Required for '.tabset'



# ##################################################

## @knitr heterogeneity_modules_expression_violin

# ##################################################

# Violinplot for module values by cluster
invisible( lapply( paste0(names(MODULES_GENES), 1), violinFeatureByCluster, seuratObject = merge_seurat_object, clustersColor = clustersColor, yLabel = "Score", addStats = FALSE, trimTitle = 1));

cat(" \n \n"); # Required for '.tabset'


