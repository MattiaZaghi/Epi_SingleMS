library(ggplot2)
library(tidyverse)
library(dplyr)

#load all the datasets to create one unique table to compare different QC Metric from the different datasets

combined.obj.ls_H3K27ac_all_mod_PT<-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_H3K27ac_all_mod_paired_tag.rds")

combined.obj.ls_H3K27me3_all_mod_PT<-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_H3K27me3_all_mod_paired_tag.rds")

combined.obj.ls_all_mod_embryo_old <-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_all_amb.rds")

combined.obj.ls_all_mod_multiome<-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_all_mod_embryo_10x.rds")

combined.obj_RNA_embryo<-list()
combined.obj_RNA_embryo[["RNA_AAAAGGGG"]]<-readRDS("/proj/user/mattia/Embryo_2/combined.obj_RNA_embryos.rds")

combined.obj_RNA_PBMCs<-list()
combined.obj_RNA_PBMCs[["RNA_AAAAGGGG"]]<-readRDS("/proj/user/mattia/Embryo_2/combined.obj_RNA_PBMCs.rds")


QC_RNA <- data.frame(
  sample = c(
    paste0(combined.obj.ls_all_mod_embryo_old[["RNA_AAAAGGGG"]]@meta.data[["modality"]], "_nanoCTAR_old"),
    paste0(combined.obj_RNA_embryo[["RNA_AAAAGGGG"]]@meta.data[["modality"]], "_nanoCTAR_new"),
    paste0(combined.obj_RNA_PBMCs[["RNA_AAAAGGGG"]]@meta.data[["modality"]], "_PBMCs_nanoCTAR_new"),
    paste0(combined.obj.ls_H3K27ac_all_mod_PT[["RNA_AAAAGGGG_H3K27ac"]]@meta.data[["modality"]], "_H3K27ac_Droplet_Pair-Tag"),
    paste0(combined.obj.ls_H3K27me3_all_mod_PT[["RNA_AAAAGGGG_H3K27me3"]]@meta.data[["modality"]], "_H3K27me3_Droplet_Pair-Tag"),
    paste0(combined.obj.ls_all_mod_multiome[["RNA_AAAAGGGG"]]@meta.data[["modality"]], "_Embryo_10x_multiome")
  ),
  UMIs = c(
    combined.obj.ls_all_mod_embryo_old[["RNA_AAAAGGGG"]]@meta.data[["nCount_RNA"]],
    combined.obj_RNA_embryo[["RNA_AAAAGGGG"]]@meta.data[["nCount_RNA"]],
    combined.obj_RNA_PBMCs[["RNA_AAAAGGGG"]]@meta.data[["nCount_RNA"]],
    combined.obj.ls_H3K27ac_all_mod_PT[["RNA_AAAAGGGG_H3K27ac"]]@meta.data[["nCount_RNA"]],
    combined.obj.ls_H3K27me3_all_mod_PT[["RNA_AAAAGGGG_H3K27me3"]]@meta.data[["nCount_RNA"]],
    combined.obj.ls_all_mod_multiome[["RNA_AAAAGGGG"]]@meta.data[["nCount_RNA"]]
  ),
  Genes = c(
    combined.obj.ls_all_mod_embryo_old[["RNA_AAAAGGGG"]]@meta.data[["nFeature_RNA"]],
    combined.obj_RNA_embryo[["RNA_AAAAGGGG"]]@meta.data[["nFeature_RNA"]],
    combined.obj_RNA_PBMCs[["RNA_AAAAGGGG"]]@meta.data[["nFeature_RNA"]],
    combined.obj.ls_H3K27ac_all_mod_PT[["RNA_AAAAGGGG_H3K27ac"]]@meta.data[["nFeature_RNA"]],
    combined.obj.ls_H3K27me3_all_mod_PT[["RNA_AAAAGGGG_H3K27me3"]]@meta.data[["nFeature_RNA"]],
    combined.obj.ls_all_mod_multiome[["RNA_AAAAGGGG"]]@meta.data[["nFeature_RNA"]]
  ),
  Mito_genes = c(
    combined.obj.ls_all_mod_embryo_old[["RNA_AAAAGGGG"]]@meta.data[["percent.mt"]],
    combined.obj_RNA_embryo[["RNA_AAAAGGGG"]]@meta.data[["percent.mt"]],
    combined.obj_RNA_PBMCs[["RNA_AAAAGGGG"]]@meta.data[["percent.mt"]],
    combined.obj.ls_H3K27ac_all_mod_PT[["RNA_AAAAGGGG_H3K27ac"]]@meta.data[["percent.mt"]],
    combined.obj.ls_H3K27me3_all_mod_PT[["RNA_AAAAGGGG_H3K27me3"]]@meta.data[["percent.mt"]],
    combined.obj.ls_all_mod_multiome[["RNA_AAAAGGGG"]]@meta.data[["percent.mt"]]
 )
)

# Define the desired order
desired_order <- c("RNA_nanoCTAR_new","RNA_PBMCs_nanoCTAR_new","RNA_nanoCTAR_old","RNA_H3K27ac_Droplet_Pair-Tag",
                   "RNA_H3K27me3_Droplet_Pair-Tag","RNA_Embryo_10x_multiome")

# Convert the sample column to a factor with the specified order
QC_RNA$sample <- factor(QC_RNA$sample, levels = desired_order)



p1=ggplot(QC_RNA) +
  aes(x = sample, y = UMIs, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)+
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 0)), 
               vjust = -5, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#3399FF", "#00CC33","#FF9933","#FFD700","#B03060")) +
  scale_color_manual(values = c("#FF0000","#3399FF", "#00CC33","#FF9933","#FFD700","#B03060")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Number of UMIs per cell") +
  scale_x_discrete(labels = c("RNA_nanoCTAR_new","RNA_PBMCs_nanoCTAR_new","RNA_nanoCTAR_old","RNA_H3K27ac_Droplet_Pair-Tag",
                              "RNA_H3K27me3_Droplet_Pair-Tag","RNA_Embryo_10x_multiome"))+
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,5000)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p1
ggsave("/proj/user/mattia/Embryo_2/plots/UMIs_plot_all.png", plot = p1, width = 20, height = 10, dpi = 300, units = "in")

p2=ggplot(QC_RNA) +
  aes(x = sample, y = Genes, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -8, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#3399FF", "#00CC33","#FF9933","#FFD700","#B03060")) +
  scale_color_manual(values = c("#FF0000","#3399FF", "#00CC33","#FF9933","#FFD700","#B03060")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Number of Genes per cell") +
  scale_x_discrete(labels = c("RNA_nanoCTAR_new","RNA_PBMCs_nanoCTAR_new","RNA_nanoCTAR_old","RNA_H3K27ac_Droplet_Pair-Tag",
                              "RNA_H3K27me3_Droplet_Pair-Tag","RNA_Embryo_10x_multiome"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,1)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p2
ggsave("/proj/user/mattia/Embryo_2/plots/Genes_plot_all.png", plot = p2, width = 20, height = 10, dpi = 300, units = "in")

p3=ggplot(QC_RNA) +
  aes(x = sample, y = Mito_genes, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -4, size = 5, color = "black")+
  scale_fill_manual(values = c("#FF0000","#3399FF", "#00CC33","#FF9933","#FFD700","#B03060")) +
  scale_color_manual(values = c("#FF0000","#3399FF", "#00CC33","#FF9933","#FFD700","#B03060")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Number of Genes per cell") +
  scale_x_discrete(labels = c("RNA_nanoCTAR_new","RNA_PBMCs_nanoCTAR_new","RNA_nanoCTAR_old","RNA_H3K27ac_Droplet_Pair-Tag",
                              "RNA_H3K27me3_Droplet_Pair-Tag","RNA_Embryo_10x_multiome"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p3
ggsave("/proj/user/mattia/Embryo_2/plots/Mito_genes_all.png", plot = p3, width = 20, height = 10, dpi = 300, units = "in")


p4=ggplot(QC_Epi) +
  aes(x = sample, y = read_per_cell, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)  +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 0)), 
               vjust = -8, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                               "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                                "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("K reads per cell") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_new","ATAC_nanoCTAR_old","ATAC_nanoCT","ATAC_Embryo_10x_multiome","H3K27ac_nanoCTAR_new","H3K27ac_nanoCTAR_old",    
                              "H3K27ac_nanoCT","H3K27ac_Droplet_Pair-Tag","H3K27ac_scCut&Tag-pro",
                              "H3K27me3_nanoCTAR_new","H3K27me3_nanoCTAR_old","H3K27me3_nanoCT","H3K27me3_Droplet_Pair-Tag","H3K27me3_scCut&Tag-pro",
                              "H3K4me3_nanoCTAR_new","H3K4me3_scCut&Tag-pro"))+
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,1)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p4
ggsave("/proj/user/mattia/Embryo_2/Reads_plot.png", plot = p4, width = 20, height = 10, dpi = 300, units = "in")




p5=ggplot(QC_Epi) +
  aes(x = sample, y = dup_rate, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)  +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -4, size = 5, color = "black")  +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                               "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                                "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Duplication rate per cell") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_new","ATAC_nanoCTAR_old","ATAC_nanoCT","ATAC_Embryo_10x_multiome","H3K27ac_nanoCTAR_new","H3K27ac_nanoCTAR_old",    
                              "H3K27ac_nanoCT","H3K27ac_Droplet_Pair-Tag","H3K27ac_scCut&Tag-pro",
                              "H3K27me3_nanoCTAR_new","H3K27me3_nanoCTAR_old","H3K27me3_nanoCT","H3K27me3_Droplet_Pair-Tag","H3K27me3_scCut&Tag-pro",
                              "H3K4me3_nanoCTAR_new","H3K4me3_scCut&Tag-pro"))+
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,1)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p5
ggsave("/proj/user/mattia/Embryo_2/dup_rate_plot.png", plot = p5, width = 20, height = 10, dpi = 300, units = "in")






