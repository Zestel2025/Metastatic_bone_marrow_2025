################################################################################
merge_2_seurat = readRDS("../p_data/merge_2_seurat_260228.rds")
merge_2_seurat$re_annotation_3_simple = substr(merge_2_seurat$re_annotation_3,start = 1,stop = 3) |> factor()
source("../script/colors.R")
levels(merge_2_seurat$re_annotation_3)

make_ids <- function(x) sprintf("c%02d", x)

merge_2_seurat$re_annotation_3_coarse <- NA
merge_2_seurat$re_annotation_3_coarse[merge_2_seurat$re_annotation_3_simple %in% make_ids(1:21)]  <- "Myeloid"
merge_2_seurat$re_annotation_3_coarse[merge_2_seurat$re_annotation_3_simple %in% make_ids(22:28)] <- "Mural"
merge_2_seurat$re_annotation_3_coarse[merge_2_seurat$re_annotation_3_simple %in% make_ids(29:32)] <- "Mesenchymal"
merge_2_seurat$re_annotation_3_coarse[merge_2_seurat$re_annotation_3_simple %in% make_ids(33:55)] <- "T cell"
merge_2_seurat$re_annotation_3_coarse[merge_2_seurat$re_annotation_3_simple %in% make_ids(56:57)] <- "NK cell"
merge_2_seurat$re_annotation_3_coarse[merge_2_seurat$re_annotation_3_simple %in% make_ids(58:65)] <- "B cell"
merge_2_seurat$re_annotation_3_coarse[merge_2_seurat$re_annotation_3_simple %in% make_ids(66)]    <- "Epithelial"

table(merge_2_seurat$re_annotation_3_coarse)

merge_2_seurat$re_annotation_3_coarse <- factor(
  merge_2_seurat$re_annotation_3_coarse,
  levels = c("Myeloid", "Mural", "Mesenchymal", "T cell", "NK cell", "B cell", "Epithelial")
)

table(merge_2_seurat$re_annotation_3_coarse)
Idents(merge_2_seurat) = merge_2_seurat$re_annotation_3_coarse
markers_All <- FindAllMarkers(object = merge_2_seurat, 
                              slot = "data",
                              test.use="wilcox",
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              min.pct = 0.1)
markers_All
write.csv(markers_All,"../Figure/Supplementary Table 1-markers of coarse clusters.csv")
Idents(merge_2_seurat) = merge_2_seurat$re_annotation_3
markers_All <- FindAllMarkers(object = merge_2_seurat, 
                              slot = "data",
                              test.use="wilcox",
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              min.pct = 0.1)
markers_All
write.csv(markers_All,"../Figure/Supplementary Table 1-markers of fine-grained clusters.csv")
