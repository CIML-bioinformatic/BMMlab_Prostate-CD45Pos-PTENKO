# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# KEGG PATHWAYS ANALYSIS
#########################

## @knitr kegg_pathway_analysis

# Get the KEGG reference
KEGG_pathway_gene_df = getGeneKEGGLinks( species="mmu")
KEGG_pathway_gene_df$Symbol = AnnotationDbi::mapIds( org.Mm.eg.db, KEGG_pathway_gene_df$GeneID, column="SYMBOL", keytype="ENTREZID")

# For each KEGG pathway to analyse perform a AddModuleScore analysis with plits and statistics
for( kegg_pathway_index in 1: nrow( KEGG_PATHWAY_DF)){
  
  kegg_pathway = KEGG_PATHWAY_DF[ kegg_pathway_index, "Pathway_name"]
  
  cat(" \n \n"); # Required for '.tabset'
  cat("### Pathway", kegg_pathway)
  cat(" \n \n"); # Required for '.tabset'
  
  # Get the KEGG ID of the pathway
  kegg_pathway_id = KEGG_PATHWAY_DF[ kegg_pathway_index, "Pathway_id"]
  
  if( ! kegg_pathway_id %in% unique( KEGG_pathway_gene_df$PathwayID)){
    cat("<BR>This pathway was not found in KEGG")
    next
  }
  
  # Get the genes symbols for the gene EntrezID
  gene_symbols = na.exclude( KEGG_pathway_gene_df[ which( KEGG_pathway_gene_df$PathwayID == paste0( "", kegg_pathway_id)), "Symbol"])
  gene_symbols = intersect( gene_symbols, rownames( sc10x.seurat))

  # Perform the AddModuleScore on the list of pathway genes
  sc10x.seurat = AddModuleScore( sc10x.seurat, features = list( gene_symbols), name = kegg_pathway)
  pathway_score_name = paste0( kegg_pathway, "1")
  
  print( FeaturePlot( sc10x.seurat, features = pathway_score_name))
  
  # Plot the modules scores by cluster and condition
  pathway_df = sc10x.seurat@meta.data[ , c( "seurat_clusters", pathway_score_name)]
  print(
    ggplot( pathway_df, 
            aes_string( x="seurat_clusters", y=pathway_score_name, fill = "seurat_clusters")) +
      stat_summary( fun = mean, geom = "point", color = "black", position = position_dodge( width = 0.9)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", size = 0.5, width = 0.3, position = position_dodge( width = 0.9)) +
      geom_hline( yintercept=0, linetype="dashed", color = "black", size = 1) + 
      theme_minimal() + ggtitle( paste( "Module scores of", kegg_pathway, "KEGG pathway (", length( gene_symbols), "genes kept)\nby cluster with mean and standard deviation"))
  ) 
}