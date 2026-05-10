source("../script/basic_plot.R")
source("../script/colors.R")
source("../script/run_barplot.R")
source("../script/propeller_function.R") # from the package speckle;
source("../script/run_stat_group.R")
# source("../script/trait_markers.R")
library(Seurat)
library(dplyr)
meta_all = readRDS("meta_all_260228.rds")
Myeloid_re_annotation_3 = readRDS("Progenitor_Myeloid_re_annotation_2.rds")
meta_Myeloid = meta_all[rownames(meta_all)%in%colnames(Myeloid_re_annotation_3),]
Myeloid_re_annotation_3@meta.data = meta_Myeloid
Myeloid_re_annotation_3$re_annotation_3_simple = substr(Myeloid_re_annotation_3$re_annotation_3,start = 1,stop = 3) %>% factor()
Idents(Myeloid_re_annotation_3)
Myeloid_re_annotation_3$re_annotation_3_simple %>% table()
Myeloid_re_annotation_3$re_annotation_3_simple = factor(Myeloid_re_annotation_3$re_annotation_3_simple)

cell_type = "Progenitor_Myeloid"
subset_name = "Myeloid_all"
seu = Myeloid_re_annotation_3
DimPlot(seu,reduction = "umap_harmony",group.by = "re_annotation_3",
        cols = col_Myeloid_HSC) + NoAxes() + ggtitle("Myeloid cell and HSC")

levels(Myeloid_re_annotation_3$re_annotation_3_simple)
all_features = c(# Progenitor
  "CD34","CRHBP","SPINK2","MSI2",#"HLF","MEIS1",
  "PCNA","TYMS","CDK6", # 23,37
  "MS4A3","PRTN3","MPO","ELANE","AZU1",# 21
  "STMN1","TOP2A","BIRC5",
  "CD33","IL17RA",
  # Neu
  "ITGAM","CEACAM8",
  "FCGR3B","MMP9","ITGAM",#"S100P",# 7,# 30
  # Mono/Mac
  "CD14","VCAN","S100A12",# 0,1,3,4,6,17,32
  "LGALS2","ISG15","SLC7A7",# 10,14,31
  "IL1B","CXCL2","NR4A2",# 2,5
  # "RNASE2","ANXA2","CES1",# 20
  #"CD68","CD163","MRC1",
  #"CCL2","CCL3", "CCL4", "CCL18",
  #"DNAJB1","HSPA1B","HSPA1A",
  "C1QB","CD68","CD163",
  "ACP5","SPP1","TREM2",# 8,9,19,22,34,38
  "FCGR3A","IFITM3","LILRB2","CX3CR1",# 11
  # Erythroid
  "GATA2","CNRIP1","HBD",#"CYTL1","PKIG",# 24
  "GATA1","CA1","REXO2",# 18
  "EPCAM","AHSP","SMIM1",#"KCNH2",# 28
  # Other
  "CD1C","FCER1A","CLEC10A",
  "CLEC9A","FLT3","IDO1",# 40
  "IRF8","SPIB","LILRB4",#pDC
  "HDC","MS4A2","TPSAB1")

pdf(file = "../Figure/Figure 2A-Dimplot_Myeloid.pdf",width = 11 ,height = 10)
DimPlot(Myeloid_re_annotation_3,reduction = "umap_harmony",group = "re_annotation_3",raster = F,
        cols = col_Myeloid_HSC ) + NoAxes() + ggtitle("Myeloid cell and HSC")
dev.off()
