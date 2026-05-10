source("../script/basic_plot.R")
source("../script/colors.R")
source("../script/run_barplot.R")
source("../script/propeller_function.R") # from the package speckle;
source("../script/run_stat_group.R")
# source("../script/trait_markers.R")
library(dplyr)
RBC_re_annotation_3 = readRDS("lineage_RBC_only.rds")
RBC_re_annotation_3$re_annotation_3_simple = substr(RBC_re_annotation_3$re_annotation_3,start = 1,stop = 3) %>% factor()
Idents(RBC_re_annotation_3) 
RBC_re_annotation_3$re_annotation_3_simple %>% table()
RBC_re_annotation_3$re_annotation_3 = factor(RBC_re_annotation_3$re_annotation_3)
RBC_re_annotation_3$re_annotation_3_simple = factor(RBC_re_annotation_3$re_annotation_3_simple)
`%nin%` <- Negate(`%in%`)
DimPlot(RBC_re_annotation_3,reduction = "umap_harmony",group.by = "harmony_cluster_3",label = T)

library(presto)
markers_All_re_annotation_1 <- FindAllMarkers(object = RBC_re_annotation_3, 
                                              slot = "data",
                                              test.use="wilcox",
                                              only.pos = TRUE,
                                              logfc.threshold = 0.25,
                                              min.pct = 0.1)

all.markers_re_annotation_1 = markers_All_re_annotation_1
all.markers_re_annotation_1[,"pct1-pct2"]=all.markers_re_annotation_1[,"pct.1"]-all.markers_re_annotation_1[,"pct.2"]
all.markers_re_annotation_1[,"abs"]=abs(all.markers_re_annotation_1[,"pct1-pct2"])
DotPlot(RBC_re_annotation_3,features = c("MME",broad),group.by = "harmony_cluster_3")+RotatedAxis()

table(RBC_re_annotation_3$cancer_type_site,RBC_re_annotation_3$re_annotation_3)
RBC_re_annotation_3_c12 = subset(RBC_re_annotation_3,subset = re_annotation_3_simple=="c15")
RBC_re_annotation_3_c13 = subset(RBC_re_annotation_3,subset = re_annotation_3_simple=="c16")
RBC_re_annotation_3_c14 = subset(RBC_re_annotation_3,subset = re_annotation_3_simple=="c17")
features = c("TFRC","HBA1","HBA2","HBB","HBD","KLF1","FTH1","FTL","TFR2")
# "FTH1": Ferritin heavy chain
# "FTL": Ferritin light chain
p_RBC_c12 = NULL
p_RBC_c13 = NULL
p_RBC_c14 = NULL
# Tumor   = "#DB3E6480",   # 80 = 大约 50% 透明
# Involved= "#E5864180",
# Distal  = "#795D9B80",
# Benign  = "#7EBAD480"
for(i in features){
  p_RBC_c12[[i]] = RidgePlot(RBC_re_annotation_3_c12, features = i,
                             ncol = 2,group.by = "cancer_site") + theme_classic()+
    scale_fill_manual(
      values = c(
        Tumor   = "#DB3E6480",   # 红
        Involved    = "#E5864180",   # 蓝
        Distal = "#795D9B80",    # 绿
        Benign = '#7EBAD480'
      ))+labs(x = paste0(i," expression"), y = "Density of cells",
              title = paste0("c12 ",i)) + NoLegend()
  df <- FetchData(
    RBC_re_annotation_3_c12,
    vars = c(i, "cancer_site", "orig.ident")
  )
  df_summary <- df %>%
    group_by(orig.ident, cancer_site) %>%
    summarise(mean_expr = mean(.data[[i]]))
  message("c12: ",i," p-value: ",kruskal.test(mean_expr ~ cancer_site, data = df_summary)$p.value)
  # 使用 pairwise.wilcox.test 进行两两比较
  pairwise_results <- pairwise.wilcox.test(
    df_summary$mean_expr,  # 要比较的值
    df_summary$cancer_site, # 分组因子
    p.adjust.method = "BH", # 多重比较校正方法
    exact = FALSE # 允许近似计算
  )
  # 查看结果
  print(pairwise_results$p.value)
  p_RBC_c13[[i]] = RidgePlot(RBC_re_annotation_3_c13, features = i,
                             ncol = 2,group.by = "cancer_site") + theme_classic()+
    scale_fill_manual(
      values = c(
        Tumor   = "#DB3E6480",   # 红
        Involved    = "#E5864180",   # 蓝
        Distal = "#795D9B80",    # 绿
        Benign = '#7EBAD480'
      ))+labs(x = paste0(i," expression"), y = "Density of cells",
              title = paste0("c13 ",i)) + NoLegend()
  df <- FetchData(
    RBC_re_annotation_3_c13,
    vars = c(i, "cancer_site", "orig.ident")
  )
  df_summary <- df %>%
    group_by(orig.ident, cancer_site) %>%
    summarise(mean_expr = mean(.data[[i]]))
  message("c13: ",i," p-value: ",kruskal.test(mean_expr ~ cancer_site, data = df_summary)$p.value)
  # 使用 pairwise.wilcox.test 进行两两比较
  pairwise_results <- pairwise.wilcox.test(
    df_summary$mean_expr,  # 要比较的值
    df_summary$cancer_site, # 分组因子
    p.adjust.method = "BH", # 多重比较校正方法
    exact = FALSE # 允许近似计算
  )
  # 查看结果
  print(pairwise_results$p.value)
  p_RBC_c14[[i]] = RidgePlot(RBC_re_annotation_3_c14, features = i,
                             ncol = 2,group.by = "cancer_site") + theme_classic()+
    scale_fill_manual(
      values = c(
        Tumor   = "#DB3E6480",   # 红
        Involved    = "#E5864180",   # 蓝
        Distal = "#795D9B80",    # 绿
        Benign = '#7EBAD480'
      ))+labs(x = paste0(i," expression"), y = "Density of cells",
              title = paste0("c14 ",i)) + NoLegend()
  df <- FetchData(
    RBC_re_annotation_3_c14,
    vars = c(i, "cancer_site", "orig.ident")
  )
  df_summary <- df %>%
    group_by(orig.ident, cancer_site) %>%
    summarise(mean_expr = mean(.data[[i]]))
  message("c14: ",i," p-value: ",kruskal.test(mean_expr ~ cancer_site, data = df_summary)$p.value)
  # 使用 pairwise.wilcox.test 进行两两比较
  pairwise_results <- pairwise.wilcox.test(
    df_summary$mean_expr,  # 要比较的值
    df_summary$cancer_site, # 分组因子
    p.adjust.method = "BH", # 多重比较校正方法
    exact = FALSE # 允许近似计算
  )
  # 查看结果
  print(pairwise_results$p.value)
}

