run_monocle3 = function(seu,remove_hsp){
  UMAP_origin = seu@reductions$umap_harmony@cell.embeddings
  # cds <- as.cell_data_set(seu)
  expression_matrix <- GetAssayData(seu, slot = "counts")
  cell_metadata <- seu@meta.data
  gene_annotation <- data.frame(
    gene_short_name = rownames(expression_matrix),
    row.names = rownames(expression_matrix)
  )
  message(Sys.time()," 正在构建cds对象……")
  cds <- new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )
  message(Sys.time()," 正在预处理...")
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- align_cds(cds, alignment_group = c("orig.ident"))
  # cds <- align_cds(cds, alignment_group = c("cancer_type"))
  cds <- reduce_dimension(cds,reduction_method = "UMAP")
  reducedDims(cds)[["UMAP"]] = UMAP_origin
  cds <- cluster_cells(cds,reduction_method = "UMAP")
  cds <- learn_graph(cds)
  # 子集提取
  message(Sys.time()," 正在提取子集...")
  cds_list <- list()
  for (g in unique(colData(cds)$cancer_type_site)) {
    message(sprintf("正在处理 %s ...", g))
    cds_subset <- tryCatch(
      {
        cds[, colData(cds)$cancer_type_site == g]
      },
      error = function(e) {
        message(sprintf("⚠️ %s 提取失败：%s", g, e$message))
        return(NULL)
      }
    )
    
    if (!is.null(cds_subset)) {
      cds_list[[g]] <- cds_subset
      message(sprintf("✅ %s 提取完成，共 %d 个细胞。", g, ncol(cds_subset)))
    }
  }
  return(list(cds = cds,
              cds_list = cds_list))
}
