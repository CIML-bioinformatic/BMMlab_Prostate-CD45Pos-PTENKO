



# FILTER DATA
# -----------

## @knitr filterData_selection

### Plot distributions of #UMIs, #Genes, %Mito, and %Ribo among cells

# Do the plots as simple images with ggplot when having a lot of points

# Generate plotly violin/jitter panels for #umis, #genes, %mitochondrial, and %ribosomal stats
ggpanel_umis  = ggplot( cbind(merge_seurat_object[["nCount_RNA"]]), aes( y = nCount_RNA)) + 
  geom_jitter( aes(x = 1.5), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
  geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
  labs(x = "# UMIs", y = "") +
  theme( axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none") +
  scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));

ggpanel_genes = ggplot( cbind(merge_seurat_object[["nFeature_RNA"]]), aes( y = nFeature_RNA)) + 
  geom_jitter( aes(x = 1.5), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
  geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
  labs(x = "# Genes", y = "") +
  theme( axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none") +
  scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));

ggpanel_mitos = if(length(merge_seurat_object$percent.mito)) ggplot( cbind(merge_seurat_object[["percent.mito"]]), aes( y = percent.mito)) + 
  geom_jitter( aes(x = 1.5), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
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


ggpanel_ribos = if(length(merge_seurat_object$percent.ribo)) ggplot( cbind(merge_seurat_object[["percent.ribo"]]), aes( y = percent.ribo)) + 
  geom_jitter( aes(x = 1.5), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
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

## @knitr filterData_summaryPlot

### Plot dispersion of excluded and non-excluded cells

# number of genes and number of UMIs
ggplot( merge_seurat_object[[]][order( merge_seurat_object[["outlier"]]),], # Plot FALSE first and TRUE after
        aes( x = nFeature_RNA, 
             y = nCount_RNA, 
             color = outlier)) + 
  geom_point( size = 0.5) +
  scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
  labs( x = "# Genes", y = "# UMIs") +
  theme( legend.position = "none")

# Mitochondrial vs ribosomal distributions
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


# Normalize and recenter the data
merge_seurat_object = NormalizeData( merge_seurat_object, normalization.method = "LogNormalize")
merge_seurat_object = ScaleData( merge_seurat_object, do.center = TRUE, do.scale = FALSE)

