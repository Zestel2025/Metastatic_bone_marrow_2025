run_nichenet = function(merge_2_seurat,
                        task = "Fibro_to_Progenitor",
                        source_id = paste0("c", 19:29),
                        target_id = paste0("c", sprintf("%02d", c(1, 3, 12, 55, 56))),
                        receiver = "c03",
                        reference = "Tumor") {
  
  cat("\n============================================================\n")
  cat("🕒 当前时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("🚀 当前任务: ", task, " | 接收细胞: ", receiver, " | 参考组: ", reference, "\n")
  cat("============================================================\n\n")
  
  # --------------------------------------------------------------------------
  cat("🕒 [", format(Sys.time(), "%H:%M:%S"), "] 1️⃣ 读取或生成 subset 数据\n")
  # --------------------------------------------------------------------------
  file_path <- paste0("./data/", task, ".rds")
  # if (file.exists(file_path)) {
  #   message("📂 读取已存在文件: ", file_path)
  #   merge_2_seurat_subset <- readRDS(file_path)
  # } else {
  #   message("⚙️ 文件不存在，开始 subset 并保存: ", file_path)
  #   merge_2_seurat_subset <- subset(
  #     merge_2_seurat,
  #     subset = (re_annotation_2_with_MP_simple %in% c(source_id, target_id)) & (orig.ident %nin% c("BMET10-Tumor"))
  #   )
  #   saveRDS(merge_2_seurat_subset, file_path)
  #   message("✅ 文件已保存: ", file_path)
  # }
  if (file.exists(file_path)) {
    message("📂 读取已存在文件: ", file_path)
    merge_2_seurat_subset <- readRDS(file_path)
    # merge_2_seurat_subset$re_annotation_2_with_MP_simple = merge_2_seurat_subset$re_annotation_2_simple
    
  } else {
    message("⚙️ 文件不存在，开始 subset 并保存: ", file_path)
    if (task == "Progenitor_to_Tc_All") {
      # 排除 BMET10-Tumor
      merge_2_seurat_subset <- subset(
        merge_2_seurat,
        subset = (re_annotation_2_with_MP_simple %in% c(source_id, target_id)) &
          !(orig.ident %in% c("BMET10-Tumor"))
      )
    } else {
      # 其他任务保留 BMET10-Tumor
      merge_2_seurat_subset <- subset(
        merge_2_seurat,
        subset = (re_annotation_2_with_MP_simple %in% c(source_id, target_id))
      )
    }
    saveRDS(merge_2_seurat_subset, file_path)
    message("✅ 文件已保存: ", file_path)
  }
  
  
  # --------------------------------------------------------------------------
  cat("🕒 [", format(Sys.time(), "%H:%M:%S"), "] 2️⃣ 预处理 Seurat 对象与基因名标准化\n")
  # --------------------------------------------------------------------------
  seuratObj = merge_2_seurat_subset
  # gene_new <- convert_alias_to_symbols(rownames(seuratObj), "human")
  # invisible(gene_new <- convert_alias_to_symbols(rownames(seuratObj), "human"))
  # dup <- duplicated(gene_new)
  # seuratObj <- seuratObj[!dup, ]
  # rownames(seuratObj) <- gene_new[!dup]
  # rownames(seuratObj) = gene_new
  
  # --------------------------------------------------------------------------
  cat("🕒 [", format(Sys.time(), "%H:%M:%S"), "] 3️⃣ 设置元信息与数据库路径\n")
  # --------------------------------------------------------------------------
  seuratObj@meta.data$celltype = seuratObj@meta.data$re_annotation_2_with_MP_simple
  seuratObj@meta.data$aggregate = seuratObj@meta.data$cancer_site
  organism <- "human"
  lr_network <- readRDS("./database/lr_network.rds")
  ligand_target_matrix <- readRDS("./database/ligand_target_matrix_nsga2r_final.rds")
  weighted_networks <- readRDS("./database/weighted_networks_nsga2r_final.rds")
  lr_network <- lr_network %>% distinct(from, to)
  # --------------------------------------------------------------------------
  cat("🕒 [", format(Sys.time(), "%H:%M:%S"), "] 4️⃣ 表达基因计算与潜在配体识别\n")
  # --------------------------------------------------------------------------
  Idents(seuratObj) = seuratObj$re_annotation_2_with_MP_simple
  expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)
  all_receptors <- unique(lr_network$to)
  expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
  potential_ligands <- lr_network %>%
    filter(to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()
  sender_celltypes <- source_id
  list_expressed_genes_sender <- sender_celltypes %>%
    unique() %>%
    lapply(get_expressed_genes, GetAssayData(seuratObj), seuratObj$celltype, 0.05)
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)
  
  # --------------------------------------------------------------------------
  cat("🕒 [", format(Sys.time(), "%H:%M:%S"), "] 5️⃣ 定义差异基因集（receiver）\n")
  # --------------------------------------------------------------------------
  if (reference == "Tumor") {
    condition_oi <- "Involved"
    condition_reference <- "Tumor"
  } else if (reference == "Involved") {
    condition_oi <- "Tumor"
    condition_reference <- "Involved"
  } else {
    stop("❌ reference must be either 'Tumor' or 'Involved'.")
  }
  seurat_obj_receiver <- subset(seuratObj, idents = receiver)
  DE_table_receiver <- FindMarkers(
    object = seurat_obj_receiver,
    ident.1 = condition_oi,
    ident.2 = condition_reference,
    group.by = "aggregate",
    min.pct = 0.05
  ) %>% rownames_to_column("gene")
  
  # geneset_oi <- DE_table_receiver %>%
  #   filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0) %>%
  #   pull(gene)
  # geneset_up <- DE_table_receiver %>%
  #   filter(p_val_adj <= 0.05 & avg_log2FC >= 0) %>% pull(gene)
  # geneset_down <- DE_table_receiver %>%
  #   filter(p_val_adj <= 0.05 & avg_log2FC <= 0) %>% pull(gene)
  # ---- Step 1: 基于padj筛选显著基因 ----
  geneset_oi <- DE_table_receiver %>%
    filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0) %>% pull(gene)
  geneset_up <- DE_table_receiver %>%
    filter(p_val_adj <= 0.05 & avg_log2FC >= 0) %>% pull(gene)
  geneset_down <- DE_table_receiver %>%
    filter(p_val_adj <= 0.05 & avg_log2FC <= 0) %>% pull(gene)
  # ---- Step 2: 如果基因数太少，退回使用未校正p值 ----
  if (length(geneset_oi) < 1000000000000000) {
    message("⚠️ 基因集太小（n=", length(geneset_oi), "），使用未校正p值进行筛选。")
    geneset_oi <- DE_table_receiver %>%
      filter(p_val <= 0.05 & abs(avg_log2FC) >= 0) %>% pull(gene)
    geneset_up <- DE_table_receiver %>%
      filter(p_val <= 0.05 & avg_log2FC >= 0) %>% pull(gene)
    geneset_down <- DE_table_receiver %>%
      filter(p_val <= 0.05 & avg_log2FC <= 0) %>% pull(gene)
  }

  geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]
  background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)]
  
  # --------------------------------------------------------------------------
  cat("🕒 [", format(Sys.time(), "%H:%M:%S"), "] 6️⃣ NicheNet ligand activity 分析\n")
  # --------------------------------------------------------------------------
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  ligand_activities <- ligand_activities %>%
    arrange(-aupr_corrected) %>%
    mutate(rank = rank(dplyr::desc(aupr_corrected)))
  
  best_upstream_ligands <- ligand_activities %>%
    top_n(300, aupr_corrected) %>%
    arrange(-aupr_corrected) %>%
    pull(test_ligand)

  
  ligand_activities_all <- ligand_activities
  best_upstream_ligands_all <- best_upstream_ligands
  
  ligand_activities <- ligand_activities %>%
    filter(test_ligand %in% potential_ligands_focused)
  best_upstream_ligands <- ligand_activities %>%
    top_n(300, aupr_corrected) %>%
    arrange(-aupr_corrected) %>%
    pull(test_ligand) %>%
    unique()
  
  ligand_aupr_matrix <- ligand_activities %>%
    filter(test_ligand %in% best_upstream_ligands) %>%
    column_to_rownames("test_ligand") %>%
    dplyr::select(aupr_corrected) %>%
    arrange(aupr_corrected)
  vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1)
  
  # --------------------------------------------------------------------------
  cat("🕒 [", format(Sys.time(), "%H:%M:%S"), "] 7️⃣ 可视化 ligand-target / receptor 网络准备\n")
  # --------------------------------------------------------------------------
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n = 3000) %>%
    bind_rows() %>%
    drop_na()
  
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0.2
  )
  
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>%
    unique() %>%
    intersect(rownames(active_ligand_target_links))
  vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])
  
  ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
    best_upstream_ligands, expressed_receptors,
    lr_network, weighted_networks$lr_sig
  )
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
    ligand_receptor_links_df,
    best_upstream_ligands,
    order_hclust = "both"
  )
  
  # --------------------------------------------------------------------------
  cat("🕒 [", format(Sys.time(), "%H:%M:%S"), "] 8️⃣ Sender-focused 差异表达（Dotplot）\n")
  # --------------------------------------------------------------------------
  seuratObj_sender = subset(seuratObj, subset = celltype %in% sender_celltypes)
  Idents(seuratObj_sender) = Idents(seuratObj_sender) %>% as.vector() %>% factor()
  celltype_order <- levels(Idents(seuratObj) %>% as.vector() %>% factor())
  DE_table_top_ligands <- lapply(
    celltype_order[celltype_order %in% sender_celltypes],
    get_lfc_celltype,
    seurat_obj = seuratObj,
    condition_colname = "aggregate",
    condition_oi = condition_oi,
    condition_reference = condition_reference,
    celltype_col = "celltype",
    min.pct = 0, logfc.threshold = 0,
    features = best_upstream_ligands
  )
  DE_table_top_ligands <- DE_table_top_ligands %>%
    purrr::reduce(., full_join) %>%
    column_to_rownames("gene")
  vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])
  
  cat("✅ 任务完成: ", task,
      " | 接收细胞: ", receiver,
      " | 参考组: ", reference,
      " | 完成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  return(list(
    vis_ligand_aupr = vis_ligand_aupr,
    seuratObj_sender = seuratObj_sender,
    best_upstream_ligands = best_upstream_ligands,
    vis_ligand_target = vis_ligand_target,
    vis_ligand_lfc = vis_ligand_lfc,
    vis_ligand_receptor_network = vis_ligand_receptor_network,
    DE_table_receiver = DE_table_receiver,
    geneset_up = geneset_up,
    geneset_down = geneset_down))
}

# run_nichenet = function(merge_2_seurat,
#                         task = "Fibro_to_Progenitor",
#                         source_id = paste0("c", 19:29),
#                         target_id = paste0("c", sprintf("%02d", c(1, 3, 12, 55, 56))),
#                         receiver = "c03",
#                         reference="Tumor"){
#   file_path <- paste0("./data/", task, ".rds")
#   if (file.exists(file_path)) {
#     # 如果文件存在，直接读取
#     message("📂 读取已存在文件: ", file_path)
#     merge_2_seurat_subset <- readRDS(file_path)
#   } else {
#     # 如果文件不存在，执行 subset 并保存
#     merge_2_seurat_subset <- subset(
#       merge_2_seurat,
#       subset = re_annotation_2_with_MP_simple %in% c(source_id, target_id))
#     saveRDS(merge_2_seurat_subset, file_path)
#     message("✅ 文件不存在，已生成并保存: ", file_path)
#   }
#   
#   seuratObj = merge_2_seurat_subset
#   
#   gene_new <- convert_alias_to_symbols(rownames(seuratObj), "human")
#   dup <- duplicated(gene_new)
#   seuratObj <- seuratObj[!dup, ]
#   rownames(seuratObj) <- gene_new[!dup]
#   rownames(seuratObj) = gene_new
#   
#   seuratObj@meta.data$celltype = seuratObj@meta.data$re_annotation_2_with_MP_simple
#   seuratObj@meta.data$aggregate = seuratObj@meta.data$cancer_site
#   organism <- "human"
#   lr_network <- readRDS("./database/lr_network.rds")
#   ligand_target_matrix <- readRDS("./database/ligand_target_matrix_nsga2r_final.rds")
#   weighted_networks <- readRDS("./database/weighted_networks_nsga2r_final.rds")
#   lr_network <- lr_network %>% distinct(from, to)
#   
#   
#   Idents(seuratObj) = seuratObj$re_annotation_2_with_MP_simple
#   expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)
#   all_receptors <- unique(lr_network$to)  
#   expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
#   potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
#   sender_celltypes <- source_id
#   # Use lapply to get the expressed genes of every sender cell type separately here
#   list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
#   expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
#   potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
#   #### 2. Define the gene set of interest ####
#   if (reference == "Tumor") {
#     condition_oi <- "Involved"
#     condition_reference <- "Tumor"
#     
#   } else if (reference == "Involved") {
#     condition_oi <- "Tumor"
#     condition_reference <- "Involved"
#     
#   } else {
#     stop("❌ reference must be either 'Tumor' or 'Involved'.")
#   }
#   seurat_obj_receiver <- subset(seuratObj, idents = receiver)
#   DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
#                                     ident.1 = condition_oi, ident.2 = condition_reference,
#                                     group.by = "aggregate",
#                                     min.pct = 0.05) %>% rownames_to_column("gene")
#   
#   geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.1) %>% pull(gene)
#   geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
#   #### 3. Define the background genes ####
#   background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
#   #### 4. Perform NicheNet ligand activity analysis #####
#   ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
#                                                  background_expressed_genes = background_expressed_genes,
#                                                  ligand_target_matrix = ligand_target_matrix,
#                                                  potential_ligands = potential_ligands)
#   ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(dplyr::desc(aupr_corrected)))
#   best_upstream_ligands <- ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
#   
#   #### 6. Sender-focused approach ####
#   ligand_activities_all <- ligand_activities 
#   best_upstream_ligands_all <- best_upstream_ligands
#   
#   ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
#   best_upstream_ligands <- ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>%
#     pull(test_ligand) %>% unique()
#   ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
#     column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
#   vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 
#   
#   
#   # Target gene plot
#   active_ligand_target_links_df <- best_upstream_ligands %>%
#     lapply(get_weighted_ligand_target_links,
#            geneset = geneset_oi,
#            ligand_target_matrix = ligand_target_matrix,
#            n = 100) %>%
#     bind_rows() %>% drop_na()
#   
#   active_ligand_target_links <- prepare_ligand_target_visualization(
#     ligand_target_df = active_ligand_target_links_df,
#     ligand_target_matrix = ligand_target_matrix,
#     cutoff = 0.2) 
#   
#   order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
#   order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
#   
#   vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
#   
#   # Receptor plot
#   ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
#     best_upstream_ligands, expressed_receptors,
#     lr_network, weighted_networks$lr_sig) 
#   
#   vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
#     ligand_receptor_links_df,
#     best_upstream_ligands,
#     order_hclust = "both") 
#   best_upstream_ligands_all %in% rownames(seuratObj) %>% table()
#   
#   # Dotplot of sender-focused approach
#   seuratObj_sender = subset(seuratObj, celltype %in% sender_celltypes)
#   Idents(seuratObj_sender) = Idents(seuratObj_sender) %>% as.vector() %>% factor()
#   
#   celltype_order <- levels(Idents(seuratObj) %>% as.vector() %>% factor()) 
#   
#   # Use this if cell type labels are the identities of your Seurat object
#   # if not: indicate the celltype_col properly
#   DE_table_top_ligands <- lapply(
#     celltype_order[celltype_order %in% sender_celltypes],
#     get_lfc_celltype, 
#     seurat_obj = seuratObj,
#     condition_colname = "aggregate",
#     condition_oi = condition_oi,
#     condition_reference = condition_reference,
#     celltype_col = "celltype",
#     min.pct = 0, logfc.threshold = 0,
#     features = best_upstream_ligands 
#   ) 
#   
#   DE_table_top_ligands <- DE_table_top_ligands %>%  purrr::reduce(., full_join) %>% 
#     column_to_rownames("gene") 
#   vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])
#   
#   return(vis_ligand_aupr = vis_ligand_aupr,
#          seuratObj_sender = seuratObj_sender,
#          vis_ligand_target = vis_ligand_target,
#          vis_ligand_receptor_network = vis_ligand_receptor_network)
# }
# 
# 

