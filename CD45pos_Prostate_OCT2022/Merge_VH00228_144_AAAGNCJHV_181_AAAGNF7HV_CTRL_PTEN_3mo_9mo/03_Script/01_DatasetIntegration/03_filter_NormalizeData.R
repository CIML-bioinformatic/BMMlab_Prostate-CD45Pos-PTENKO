
# FILTER DATA
# -----------

## @knitr filterData_selection


if( QC_EXPLORATION_MODE == TRUE){
  cat("<p style='color:red;'><b>WARNING: QC EXPLORATION MODE is active : Cells will not be filtered but only marked as different</b></p>")
}

### Identify mitocondrial genes in matrix
mito.genes = grep( pattern = "^mt-", x = rownames( LayerData(object = merge_seurat_object, layer = "counts")), value = TRUE, ignore.case = TRUE)
if(length(mito.genes)==0) {
  warning( "No mitochondrial genes could be identified in this dataset.");
} else {
  # Compute percentage of mitochondrial transcripts in each cell
  percent.mito <- PercentageFeatureSet(merge_seurat_object, features=mito.genes)
  # Add the mitocondrial gene percentage as meta information in the Seurat object
  merge_seurat_object[["percent.mito"]] <- percent.mito
}

### Identify ribosomal genes in matrix
ribo.genes = grep(pattern = "^rps|^rpl",  x = rownames(GetAssayData(object = merge_seurat_object, slot = "counts")), value = TRUE, ignore.case=TRUE)
if(length(ribo.genes)==0) {
  warning( "No ribosomal genes could be identified in this dataset.");
} else {
  # Compute percentage of ribosomal transcripts in each cell
  percent.ribo <- PercentageFeatureSet(merge_seurat_object, features=ribo.genes)
  # Add the ribosomal gene percentage as meta information in the Seurat object
  merge_seurat_object[["percent.ribo"]] <- percent.ribo
}


### Identify cells that will be rejected based on specified thresholds

# Reject cells based on UMI numbers
nUMI.drop = logical( length(Cells(merge_seurat_object)));
if( ! is.null( FILTER_UMI_MIN)){
  nUMI.drop = nUMI.drop | (merge_seurat_object[["nCount_RNA", drop=TRUE]] < FILTER_UMI_MIN);
}

if( ! is.null( FILTER_UMI_MAX)){
  nUMI.drop = nUMI.drop | (merge_seurat_object[["nCount_RNA", drop=TRUE]] > FILTER_UMI_MAX);
}
merge_seurat_object = AddMetaData( merge_seurat_object, metadata = nUMI.drop[ Cells( merge_seurat_object)], col.name = "outlier.nCount_RNA")

# Reject cells based on number of expressed genes
nGene.drop = logical( length(Cells(merge_seurat_object)));
if( ! is.null( FILTER_FEATURE_MIN)){
  nGene.drop = nGene.drop | (merge_seurat_object[["nFeature_RNA", drop=TRUE]] < FILTER_FEATURE_MIN);
}

if( ! is.null( FILTER_FEATURE_MAX)){
  nGene.drop = nGene.drop | (merge_seurat_object[["nFeature_RNA", drop=TRUE]] > FILTER_FEATURE_MAX);
}
merge_seurat_object = AddMetaData( merge_seurat_object, metadata = nGene.drop[ Cells( merge_seurat_object)], col.name = "outlier.nFeature_RNA")

# Identify cells with high percentage of mitocondrial genes
mito.drop = logical( length(Cells(merge_seurat_object)));
if( length(mito.genes) && (! is.null(FILTER_MITOPCT_MAX))){
  mito.drop = (merge_seurat_object[["percent.mito", drop=TRUE]] > FILTER_MITOPCT_MAX);
}
merge_seurat_object = AddMetaData( merge_seurat_object, metadata = mito.drop[ Cells( merge_seurat_object)], col.name = "outlier.percent.mito")

# Identify cells with low percentage of ribosomal genes
ribo.drop = logical( length(Cells(merge_seurat_object)));
if( length(ribo.genes) && (! is.null(FILTER_RIBOPCT_MIN))){
  ribo.drop = (merge_seurat_object[["percent.ribo", drop=TRUE]] < FILTER_RIBOPCT_MIN);
}
merge_seurat_object = AddMetaData( merge_seurat_object, metadata = ribo.drop[ Cells( merge_seurat_object)], col.name = "outlier.percent.ribo")


### Plot distributions of #UMIs, #Genes, %Mito, and %Ribo among cells

# Do the plots as simple images with ggplot when having a lot of points

# Generate plotly violin/jitter panels for #umis, #genes, %mitochondrial, and %ribosomal stats
ggpanel_umis  = ggplot( cbind(merge_seurat_object[["nCount_RNA"]], outliers = nUMI.drop), aes( y = nCount_RNA)) + 
  geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
  geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
  geom_hline( yintercept = c( FILTER_UMI_MIN, FILTER_UMI_MAX), 
              color = c( "blue", "red")[!sapply( list( FILTER_UMI_MIN, FILTER_UMI_MAX), is.null)], 
              alpha = 0.5,
              size = 1) +
  labs(x = "# UMIs", y = "") +
  theme( axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none") +
  scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));

ggpanel_genes = ggplot( cbind(merge_seurat_object[["nFeature_RNA"]], outliers = nGene.drop), aes( y = nFeature_RNA)) + 
  geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
  geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
  geom_hline( yintercept = c( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), 
              color = c( "blue", "red")[!sapply( list( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), is.null)], 
              alpha = 0.5,
              size = 1) +
  labs(x = "# Genes", y = "") +
  theme( axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none") +
  scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));

ggpanel_mitos = if(length(mito.genes)) ggplot( cbind(merge_seurat_object[["percent.mito"]], outliers = mito.drop), aes( y = percent.mito)) + 
  geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
  geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
  geom_hline( yintercept = FILTER_MITOPCT_MAX, 
              color = "red", 
              alpha = 0.5,
              size = 1) +
  labs(x = "% Mito", y = "") +
  theme( axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none") +
  scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444")) else NULL;


ggpanel_ribos = if(length(ribo.genes)) ggplot( cbind(merge_seurat_object[["percent.ribo"]], outliers = ribo.drop), aes( y = percent.ribo)) + 
  geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
  geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
  geom_hline( yintercept = FILTER_RIBOPCT_MIN, 
              color = "blue", 
              alpha = 0.5,
              size = 1) +
  labs(x = "% Ribo", y = "") +
  theme( axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none") +
  scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444")) else NULL;

# Use patchwork library to assemble panels
print( ggpanel_umis + ggpanel_genes + ggpanel_mitos + ggpanel_ribos + plot_layout( nrow = 1))



cat( "<br>Number of cells removed based on number of UMIs:", sum( nUMI.drop));
cat( "<br>Number of cells removed based on number of genes:", sum( nGene.drop));
if(exists( "mito.drop")) cat( "<br>Number of cells removed based on high percentage of mitochondrial transcripts:", sum( mito.drop));
if(exists( "ribo.drop")) cat( "<br>Number of cells removed based on low percentage of ribosomal transcripts:", sum( ribo.drop));
cat( "\n<br>\n");

# Identify cells to exclude as union of cells with low nb UMI, low nb expressed genes, high pct mito genes, low pct ribo genes
merge_seurat_object[["outlier"]] = nUMI.drop  | 
  nGene.drop | 
  ( if(exists( "mito.drop")) mito.drop else FALSE ) | 
  ( if(exists( "ribo.drop")) ribo.drop else FALSE );

cat("<br><br>Removed cells after filters:", sum( unlist(merge_seurat_object[["outlier"]] )));
cat("<br>Remaining cells after filters:", sum( ! unlist(merge_seurat_object[["outlier"]] )));
cat("\n<br>\n");

### Record which cells got rejected

# Export the excluded cells to file
write.table( data.frame( cellid = names(which(nUMI.drop))), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, 
                              paste0( outputFilesPrefix, "excluded_cells_nbUMI.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");

write.table( data.frame( cellid = names(which(nGene.drop))), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_nbGene.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");

if(exists( "mito.drop")){
  write.table( data.frame( cellid = names(which(mito.drop))), 
               file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctMito.txt")), 
               quote = FALSE, 
               row.names = FALSE, 
               col.names = TRUE, 
               sep="\t");
}

if(exists( "ribo.drop")){
  write.table( data.frame( cellid = names(which(ribo.drop))), 
               file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctRibo.txt")), 
               quote = FALSE, 
               row.names = FALSE, 
               col.names = TRUE, 
               sep="\t");
}




## @knitr filterData_summaryPlot

### Plot dispersion of excluded and non-excluded cells

# number of genes and number of UMIs
ggplot( merge_seurat_object[[]][order( merge_seurat_object[["outlier"]]),], # Plot FALSE first and TRUE after
        aes( x = nFeature_RNA, 
             y = nCount_RNA, 
             color = outlier)) + 
  geom_point( size = 0.5) +
  geom_vline( xintercept = c( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), 
              linetype = 2, 
              color = c( "blue", "red")[!sapply( list( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), is.null)], 
              alpha = 0.5) +
  geom_hline( yintercept = c( FILTER_UMI_MIN, FILTER_UMI_MAX), 
              linetype = 2, 
              color = c( "blue", "red")[!sapply( list( FILTER_UMI_MIN, FILTER_UMI_MAX), is.null)], 
              alpha = 0.5) +
  scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
  labs( x = "# Genes", y = "# UMIs") +
  theme( legend.position = "none")

# Mitochondrial vs ribosomal distributions
if(exists( "mito.drop") && exists( "ribo.drop")){
  ggplot( merge_seurat_object[[]][order( merge_seurat_object[["outlier"]]),], # Plot FALSE first and TRUE after
          aes( x = percent.ribo, 
               y = percent.mito, 
               color = outlier)) + 
    geom_point( size = 0.5) +
    geom_vline( xintercept = FILTER_RIBOPCT_MIN, linetype = 2, color = "blue", alpha = 0.5) +
    geom_hline( yintercept = FILTER_MITOPCT_MAX, linetype = 2, color = "red", alpha = 0.5) +
    scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
    labs( x = "% Ribosomal genes", y = "% Mitochondrial genes") +
    theme( legend.position = "none")
}


# FILTER DATA
# --------------

## @knitr filterData_filterObject

# Filter the excluded cells in the Seurat object
if( QC_EXPLORATION_MODE == FALSE){
  merge_seurat_object = merge_seurat_object[ , ! merge_seurat_object[[ "outlier", drop=TRUE ]] ];
}
# Save the list of remaining cells after selection during loading, and #Genes, #UMIs, pctMito, pctRibo
write.table( data.frame( cellid = Cells(merge_seurat_object)), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "selected_cells.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");


# DETECT DOUBLET
# --------------

## @knitr detect_doublet

### Read in data as an sce object ###
sce <- SingleCellExperiment( list( counts= merge_seurat_object@assays$RNA@counts))

## Calculate doublet ratio ###
doublet_ratio <- ncol(sce)/1000*0.008

### Calculate Singlets and Doublets ###
sce <- scDblFinder(sce, dbr=doublet_ratio)

cell_class = sce$scDblFinder.class
names( cell_class) = colnames( sce)

merge_seurat_object = AddMetaData( merge_seurat_object, metadata = cell_class[ Cells( merge_seurat_object)], col.name = "scDblFinder")

# Show the table of doublet assignation against outliers
if( QC_EXPLORATION_MODE == TRUE){ 
  cat("<BR>")
  table( merge_seurat_object$outlier, merge_seurat_object$scDblFinder) %>% kable() %>% kable_styling( full_width = FALSE)
  cat("<BR>")
}

# Show the table of doublet assignation against sample of origin
cat("<BR>")
table( merge_seurat_object$scDblFinder, merge_seurat_object$orig.ident) %>% kable() %>% kable_styling( full_width = FALSE)
cat("<BR>")


# NORMALIZE DATA
# --------------

## @knitr normalizeData

merge_seurat_object = NormalizeData( object = merge_seurat_object,
                       normalization.method = DATA_NORM_METHOD,
                       scale.factor = DATA_NORM_SCALEFACTOR,
                       verbose = .VERBOSE)

merge_seurat_object = ScaleData( object    = merge_seurat_object,
                                 do.center = DATA_CENTER,
                                 do.scale  = DATA_SCALE,
                                 vars.to.regress = DATA_VARS_REGRESS,
                                 verbose = .VERBOSE)

