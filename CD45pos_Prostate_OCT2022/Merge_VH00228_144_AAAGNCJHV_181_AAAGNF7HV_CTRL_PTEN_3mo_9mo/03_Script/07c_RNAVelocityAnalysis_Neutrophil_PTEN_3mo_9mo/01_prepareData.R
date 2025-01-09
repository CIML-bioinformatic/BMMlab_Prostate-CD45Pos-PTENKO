# ################################################
# This script aims to read and filter data
# ################################################

# This script has been inspired by the R Velocyto tutorial page:
# http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html

## @knitr prepare_data

# READ THE DATA
# ...............

## @knitr prepare_data

# Load the Seurat object from heterogenity analysis
sc10x.seurat = readRDS( SEURAT_OBJECT_HETEROGENEITY_ANALYSIS)

# Select the cells on condition, date and clusters
selected_cells_condition = Cells( sc10x.seurat)[ sc10x.seurat[[ "orig.condition"]] == SELECTED_CONDITION]
selected_cells_cluster = Cells( sc10x.seurat)[ which( sc10x.seurat[[ "seurat_clusters"]]$seurat_clusters %in% CLUSTER_GROUP_LIST[[ CLUSTER_GROUP]] )]
selected_cells = intersect( selected_cells_condition, selected_cells_cluster)

cat("<BR><b>WARNING : Selecting cells from condition", SELECTED_CONDITION, "in clusters", paste( CLUSTER_GROUP_LIST[[ CLUSTER_GROUP]], collapse = ',') , "</b>")
cat("<BR>Number of cells of condition", SELECTED_CONDITION,"=", length( selected_cells_condition))
cat("<BR>Number of cells of clusters", paste( CLUSTER_GROUP_LIST[[ CLUSTER_GROUP]], collapse = ','),"=", length( selected_cells_cluster))
cat("<BR>Number of selected cells =", length( selected_cells))
cat("<BR>Number of excluded cells =", length( Cells( sc10x.seurat)) - length( selected_cells))

# Subset the seurat object with the selected the cells and relevel the seurat clusters factor
sc10x.seurat = subset( sc10x.seurat, cells = selected_cells)
sc10x.seurat$seurat_clusters = factor( sc10x.seurat$seurat_clusters, levels=sort( unique( sc10x.seurat$seurat_clusters)))

# Show the distribution of selected cell in clusters
table( sc10x.seurat$seurat_clusters, sc10x.seurat$orig.condition) %>% 
            kable( caption = "Number of cells per cluster", col.names = c( "Cluster", "Number of cells")) %>% 
            kable_styling( full_width = FALSE)

# Plot the UMAP split by condition and colored by clusters
DimPlot(sc10x.seurat, group.by = "seurat_clusters", split.by = "orig.condition")

# Show the dispersion of cell type prediction (L1) to check for the dataset composition
table( sc10x.seurat$predicted.celltype.l1) %>% kable( caption = "Dispersion of predicted L1 cell types") %>% kable_styling( full_width = FALSE)

# Get the UMAP, clusters and cells by date from the Seurat object
umap_embedding = as.data.frame( sc10x.seurat@reductions$umap@cell.embeddings)
umap_embedding$cluster = sc10x.seurat[[ "seurat_clusters"]][ row.names( umap_embedding), "seurat_clusters"]
umap_embedding$date = sc10x.seurat[[ "orig.date"]][ row.names( umap_embedding), "orig.date"]


# Read the loom file from velocyto pre-processing
input_file = file.path( VELOCYTO_LOOM_FILE)
ldat <- read.loom.matrices( input_file)
cat("<BR>Analyzed data are located on:", input_file)

# Get expression matrix
emat <- ldat$spliced

# Convert the names of the cells (they are different in the Seurat Object and in the loom object)
# Add the prefix of the cell in the loom to the cell name in the UMAP embedding obtained from seurat object
# Convert the "-1" (seurat style) at the end of cell name to "x" (loom style)
emat_cell_name = colnames( emat)[1]
emat_prefix = str_match(emat_cell_name, "(.*):(.*)")[2]
row.names( umap_embedding) = gsub( "-1", "x", row.names( umap_embedding))
row.names( umap_embedding) = paste0( emat_prefix, ":" ,row.names( umap_embedding))

# Filter some very low values
emat = emat[ , row.names( umap_embedding)]
#emat <- emat[,colSums(emat)>=1e3]


# Pagoda2 pre-processing
# .......................

## @knitr build_data
cat("<BR>Building the required data for velocity analysis")

# Create the Pagoda2 object
rownames( emat) <- make.unique( rownames( emat))
cat("<BR>Number of genes in expression matrix:", length( rownames( emat)))

r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)

# Adjust the variances
r$adjustVariance(plot=T,do.par=T,gam.k=10)

# Run basic analysis steps to generate cell embedding and clustering, visualize:

# - Run the PCA
set.seed( 123)
cat("<HR>")
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
cat("<HR>")

# - Cluster the cells
r$makeKnnGraph(k=20,type='PCA',center=T,distance='cosine');
r$getKnnClusters(method=multilevel.community, type='PCA', name='multilevel')
#r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=FALSE)

seurat_cluster_set = umap_embedding[ names( r$clusters[[ "PCA"]]$multilevel), "cluster"]
names( seurat_cluster_set) = names( r$clusters[[ "PCA"]]$multilevel)

r$clusters[[ "PCA"]]$seurat_clusters = seurat_cluster_set

# - Look at the distribution of cells along cluster
cluster_df = data.frame( table( r$clusters[[ "PCA"]]$seurat_clusters))
names( cluster_df) = c( "Cluster ID", "Cell Number")
cluster_df %>% kbl( caption="Distribution of cells along clusters", rownames=FALSE) %>%kable_styling( full_width = FALSE)

# Compute differentially expressed genes
deg_gene_list = r$getDifferentialGenes(type='PCA', verbose=FALSE, clusterType='seurat_clusters')
de = data.frame()
for( cluster_name in names( r$diffgenes[[ "PCA"]]$seurat_clusters)){
  current_df = head( r$diffgenes[[ "PCA"]]$seurat_clusters[[ cluster_name]], 10)
  current_df$cluster = rep( paste0( "cluster", cluster_name), nrow( current_df))
  de = rbind( de, current_df)
}
cat("<HR>")
#r$plotGeneHeatmap(genes = rownames( de), groups = r$clusters$PCA[[1]])
cat("<HR>")
de %>% kbl( caption = "Differentially expressed genes by cluster") %>% kable_styling( full_width = FALSE)
cat("<HR>")

# Plot embedding, labeling clusters
# par( mfrow = c(1,1))
# cat("<HR>")
# r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
# cat("<HR>")