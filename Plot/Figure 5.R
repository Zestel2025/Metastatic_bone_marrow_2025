p1 = FeaturePlot(CD4_re_annotation_2, features = "TCF7", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 = FeaturePlot(CD4_re_annotation_2, features = "CCR7", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p3 = FeaturePlot(CD4_re_annotation_2, features = "TIGIT", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p4 = FeaturePlot(CD4_re_annotation_2, features = "LPAR6", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("") + scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p5 = FeaturePlot(CD4_re_annotation_2, features = "KLRB1", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+ scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p6 = FeaturePlot(CD4_re_annotation_2, features = "CD69", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p7 = FeaturePlot(CD4_re_annotation_2, features = "GZMA", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p8 = FeaturePlot(CD4_re_annotation_2, features = "ISG15", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p9 = FeaturePlot(CD4_re_annotation_2, features = "FOXP3", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p10 = FeaturePlot(CD4_re_annotation_2, features = "NKG7", max.cutoff = 'q95',reduction = "umap_harmony",order = T) + 
  NoLegend() + NoAxes()+ ggtitle("")+scale_color_viridis() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
