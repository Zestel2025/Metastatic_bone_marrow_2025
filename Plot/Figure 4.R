source("../script/run_barplot.R")
cell_type = "re_annotation_2_simple"
all_proportion = run_barplot(Bc_re_annotation_2,
                             group = c("Benign","Distal","Involved","Tumor"),
                             cell_type = cell_type)
all_proportion$celltype = factor(all_proportion$celltype,levels = levels(Bc_re_annotation_2@meta.data[,cell_type]))
pdf("../Figure/Figure S4-Barplot-Bc.pdf",width = 9, height = 6)
ggplot(all_proportion, aes(x=celltype, y=mean, fill=clade)) + 
  geom_bar(stat="identity", color="black", position = position_dodge()) + 
  theme_classic() + xlab("") + ylab("") + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
  theme(legend.position = 'none', 
        # axis.line = element_line(size=1.25),
        # axis.ticks = element_line(size=1.25),
        # axis.ticks.length=unit(.2, "cm"),
        axis.text.x = element_text(# family = "Arial", 
          color = "black", size = 15), 
        axis.text.y = element_text(# family = "Arial", 
          color = "black", size = 15)
  ) + 
  scale_fill_manual(values=c('#7EBAD4','#795D9B','#E58641','#DB3E64')) +RotatedAxis()
dev.off()
