# ##################################################
# This script compare the distribution of cells
# across the co-variables
# ##################################################



# .....................................................
## @knitr sample_vs_cluster_distribution
# .....................................................

# Look at the number of cells per samples and cluster
Idents( merge_seurat_object) = "seurat_clusters"
cells_clusterid = Idents( merge_seurat_object)
Idents( merge_seurat_object) = "orig.ident"
cells_sample = Idents( merge_seurat_object)

# Look at the distribution of sample across clusters
sample_vs_cluster_table = table( cells_clusterid, cells_sample)
sample_vs_cluster_chisq_test = chisq.test( sample_vs_cluster_table)

sample_vs_cluster_df = as.data.frame.table( sample_vs_cluster_table)
ggplot( sample_vs_cluster_df, aes(fill=cells_sample, y=Freq, x=cells_clusterid)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual( values = SAMPLE_COLOR)

sample_vs_cluster_table %>% kable( caption = "Sample of origin versus cluster") %>% kable_styling( full_width = FALSE)
corrplot::corrplot( sample_vs_cluster_chisq_test$residuals, is.cor = FALSE)

# .....................................................
## @knitr condition_vs_cluster_distribution
# .....................................................

# Look at the number of cells per condition and cluster
Idents( merge_seurat_object) = "seurat_clusters"
cells_clusterid = Idents( merge_seurat_object)
Idents( merge_seurat_object) = "orig.condition"
cells_condition = Idents( merge_seurat_object)

# Look at the distribution of condition across clusters
condition_vs_cluster_table = table( cells_clusterid, cells_condition)
condition_vs_cluster_chisq_test = chisq.test( condition_vs_cluster_table)

condition_vs_cluster_df = as.data.frame.table( condition_vs_cluster_table)
ggplot( condition_vs_cluster_df, aes(fill=cells_condition, y=Freq, x=cells_clusterid)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual( values = CONDITION_COLOR)

condition_vs_cluster_table %>% kable( caption = "condition of origin versus cluster") %>% kable_styling( full_width = FALSE)
corrplot::corrplot( condition_vs_cluster_chisq_test$residuals, is.cor = FALSE)

# .....................................................
## @knitr cluster_vs_condition_distribution
# .....................................................

# Look at the number of cells per condition and cluster
Idents( merge_seurat_object) = "seurat_clusters"
cells_clusterid = Idents( merge_seurat_object)
Idents( merge_seurat_object) = "orig.condition"
cells_condition = Idents( merge_seurat_object)

# Look at the distribution of condition across clusters
condition_vs_cluster_table = table( cells_clusterid, cells_condition)
condition_vs_cluster_chisq_test = chisq.test( condition_vs_cluster_table)

condition_vs_cluster_df = as.data.frame.table( condition_vs_cluster_table)
ggplot( condition_vs_cluster_df, aes(fill=cells_clusterid, y=Freq, x=cells_condition)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual( values = CLUSTER_COLOR_PANEL[ 1 : length( levels( cells_clusterid))])

condition_vs_cluster_table %>% kable( caption = "condition of origin versus cluster") %>% kable_styling( full_width = FALSE)
corrplot::corrplot( condition_vs_cluster_chisq_test$residuals, is.cor = FALSE)

# .....................................................
## @knitr condition_vs_celltypeL1_distribution
# .....................................................

# Look at the number of cells per condition and cell type
Idents( merge_seurat_object) = "predicted.celltype.l1"
cells_celltype = Idents( merge_seurat_object)
Idents( merge_seurat_object) = "orig.condition"
cells_condition = Idents( merge_seurat_object)

# Look at the distribution of condition across celltypes
condition_vs_celltype_table = table( cells_celltype, cells_condition)
condition_vs_celltype_chisq_test = chisq.test( condition_vs_celltype_table)

condition_vs_celltype_df = as.data.frame.table( condition_vs_celltype_table)
ggplot( condition_vs_celltype_df, aes(fill=cells_condition, y=Freq, x=cells_celltype)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual( values = CONDITION_COLOR)

condition_vs_celltype_table %>% kable( caption = "condition of origin versus cell type") %>% kable_styling( full_width = FALSE)
corrplot::corrplot( condition_vs_celltype_chisq_test$residuals, is.cor = FALSE)

# .....................................................
## @knitr clusters_vs_celltypeL1_distribution
# .....................................................

# Look at the number of cells per cluster and cell type
Idents( merge_seurat_object) = "predicted.celltype.l1"
cells_celltype = Idents( merge_seurat_object)
Idents( merge_seurat_object) = "seurat_clusters"
cells_cluster = Idents( merge_seurat_object)

# Look at the distribution of cluster across celltypes
cluster_vs_celltype_table = table( cells_cluster, cells_celltype)
cluster_vs_celltype_chisq_test = chisq.test( cluster_vs_celltype_table)

cluster_vs_celltype_df = as.data.frame.table( cluster_vs_celltype_table)
ggplot( cluster_vs_celltype_df, aes(fill=cells_celltype, y=Freq, x=cells_cluster)) + 
  geom_bar(position="fill", stat="identity")

cluster_vs_celltype_table %>% kable( caption = "cluster versus cell type") %>% kable_styling( full_width = FALSE)
corrplot::corrplot( cluster_vs_celltype_chisq_test$residuals, is.cor = FALSE)

