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
## @knitr date_vs_cluster_distribution
# .....................................................

# Look at the number of cells per date and cluster
Idents( merge_seurat_object) = "seurat_clusters"
cells_clusterid = Idents( merge_seurat_object)
Idents( merge_seurat_object) = "orig.date"
cells_date = Idents( merge_seurat_object)

# Look at the distribution of date across clusters
date_vs_cluster_table = table( cells_clusterid, cells_date)
date_vs_cluster_chisq_test = chisq.test( date_vs_cluster_table)

date_vs_cluster_df = as.data.frame.table( date_vs_cluster_table)
ggplot( date_vs_cluster_df, aes(fill=cells_date, y=Freq, x=cells_clusterid)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual( values = DATE_COLOR)

date_vs_cluster_table %>% kable( caption = "date of origin versus cluster") %>% kable_styling( full_width = FALSE)
corrplot::corrplot( date_vs_cluster_chisq_test$residuals, is.cor = FALSE)

# .....................................................
## @knitr cluster_vs_date_distribution
# .....................................................

# Look at the number of cells per date and cluster
Idents( merge_seurat_object) = "seurat_clusters"
cells_clusterid = Idents( merge_seurat_object)
Idents( merge_seurat_object) = "orig.date"
cells_date = Idents( merge_seurat_object)

# Look at the distribution of date across clusters
date_vs_cluster_table = table( cells_clusterid, cells_date)
date_vs_cluster_chisq_test = chisq.test( date_vs_cluster_table)

date_vs_cluster_df = as.data.frame.table( date_vs_cluster_table)
ggplot( date_vs_cluster_df, aes(fill=cells_clusterid, y=Freq, x=cells_date)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual( values = CLUSTER_COLOR_PANEL[ 1 : length( levels( cells_clusterid))])

date_vs_cluster_table %>% kable( caption = "date of origin versus cluster") %>% kable_styling( full_width = FALSE)
corrplot::corrplot( date_vs_cluster_chisq_test$residuals, is.cor = FALSE)

# .....................................................
## @knitr date_vs_celltypeL1_distribution
# .....................................................

# Look at the number of cells per date and cell type
Idents( merge_seurat_object) = "predicted.celltype.l1"
cells_celltype = Idents( merge_seurat_object)
Idents( merge_seurat_object) = "orig.date"
cells_date = Idents( merge_seurat_object)

# Look at the distribution of date across celltypes
date_vs_celltype_table = table( cells_celltype, cells_date)
date_vs_celltype_chisq_test = chisq.test( date_vs_celltype_table)

date_vs_celltype_df = as.data.frame.table( date_vs_celltype_table)
ggplot( date_vs_celltype_df, aes(fill=cells_date, y=Freq, x=cells_celltype)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual( values = DATE_COLOR)

date_vs_celltype_table %>% kable( caption = "date of origin versus cell type") %>% kable_styling( full_width = FALSE)
corrplot::corrplot( date_vs_celltype_chisq_test$residuals, is.cor = FALSE)

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

