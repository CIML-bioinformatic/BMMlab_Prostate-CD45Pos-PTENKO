
SCORE_MODE == "EXPRESSION"

if( SCORE_MODE == "EXPRESSION"){
goterm_score_df = data.frame()
for( goterm in selected_bp_go_onto_df$id){
  
  gene_list = str_split( selected_bp_go_onto_df[ goterm, "gene.list"], pattern = ";")[[1]]
  gene_list = intersect( gene_list, rownames( merge_seurat_object))
  
  if( length( gene_list) > 10){
    goterm_score_df = rbind( goterm_score_df, t( data.frame( colSums( merge_seurat_object@assays$RNA@counts[ gene_list,]))))
  }
}
row.names( goterm_score_df) = gsub( ":", "-", selected_bp_go_onto_df$id)

goterm_seurat = CreateSeuratObject( goterm_score_df, project = "GoTermProject", assay = "RNA",
                                    min.cells = 0, min.features = 0, names.field = 1,
                                    names.delim = "_", meta.data = NULL)
}else{
  features_list = list()
  for( goterm in selected_bp_go_onto_df$id){
    
    gene_list = str_split( selected_bp_go_onto_df[ goterm, "gene.list"], pattern = ";")[[1]]
    gene_list = intersect( gene_list, rownames( merge_seurat_object))
    
    if( length( gene_list) > 10){
      features_list[[ gsub( ":", "-", goterm)]] = gene_list
    }
  }
  featured_merge_seurat_object = AddModuleScore( object = merge_seurat_object, features = features_list)
  
  goterm_score_df = t( featured_merge_seurat_object@meta.data[ , paste0( "Cluster", seq( 1, 2449,1))])
  row.names( goterm_score_df) = names( features_list)
  
  goterm_seurat = CreateSeuratObject( goterm_score_df, project = "GoTermProject", assay = "RNA",
                                      min.cells = 0, min.features = 0, names.field = 1,
                                      names.delim = "_", meta.data = NULL)
}

# Normalize and recenter the data
goterm_seurat = NormalizeData( goterm_seurat, normalization.method = "LogNormalize")
goterm_seurat@assays$RNA@layers$data = goterm_seurat@assays$RNA@layers$counts
goterm_seurat = ScaleData( goterm_seurat, do.center = TRUE, do.scale = TRUE, assay = "RNA")

goterm_seurat <- RunPCA( object = goterm_seurat,
                         features = Features( goterm_seurat),
                         npcs     = 30,
                         verbose  = .VERBOSE)

# Run the UMAP with the optimized parameters
goterm_seurat = RunUMAP( goterm_seurat, 
                         dims = 1:30, 
                         n.neighbors = 30, 
                         min.dist = 0.3)

Idents( merge_seurat_object) = "seurat_clusters"
goterm_seurat = AddMetaData( object = goterm_seurat, metadata = Idents( merge_seurat_object), col.name = "RNA_seurat_clusters")

Idents( merge_seurat_object) = "orig.condition"
goterm_seurat = AddMetaData( object = goterm_seurat, metadata = Idents( merge_seurat_object), col.name = "orig.condition")

# Plot the cells clusters issued from RNA analysis
DimPlot( goterm_seurat, reduction = "umap", group.by = "RNA_seurat_clusters") + 
  ggtitle( "Map of cells by clusters from RNA study")
DimPlot( goterm_seurat, reduction = "umap", group.by = "orig.condition", split.by = "orig.condition", ncol = 2) + 
  ggtitle( "Map of cells by condition")
DimPlot( merge_seurat_object, reduction = "umap", group.by = "orig.condition", split.by = "orig.condition", ncol = 2) + 
  ggtitle( "Map of cells by condition")

goterm_seurat <- FindNeighbors(object    = goterm_seurat,
                                     k.param   = 30, 
                                     reduction = "pca",
                                     dims      = 1:30,
                                     verbose   = .VERBOSE);

goterm_seurat = FindClusters( object             = goterm_seurat,
                              resolution         = 0.8,
                              algorithm          = 1,
                              temp.file.location = "/tmp/",
                              verbose            = FALSE)

DimPlot( goterm_seurat, reduction = "umap",
         order = TRUE, 
         label = TRUE, 
         label.size = 6)  +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none")

DimPlot( goterm_seurat, reduction = "umap", split.by = "orig.condition", ncol = 2,
         order = TRUE, 
         label = TRUE, 
         label.size = 6)  +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none")

# Identify marker genes
Idents( goterm_seurat) = "seurat_clusters"
markers = FindAllMarkers( object          = goterm_seurat,
                          test.use        = "roc",
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          return.thresh   = 0.05,
                          random.seed     = SEED,
                          verbose         = .VERBOSE)

markers$GoTermName = selected_bp_go_onto_df[ gsub( "-", ":", markers$gene), "name"]

FeaturePlot( goterm_seurat, features = c( "GO-0048701"), split.by = "orig.condition", ncol = 2 )

VlnPlot(object = goterm_seurat, features = c( "GO-0048701"), split.by = "orig.condition")  

markers[ which( markers$cluster == "4"), ]
