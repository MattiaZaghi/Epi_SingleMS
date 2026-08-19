library(ggplot2)
library(tidyverse)
library(dplyr)
library(Seurat)
library(Signac)

obj.ls.qc_nanoCTAR<-readRDS("/proj/user/mattia/Embryo_2/obj.ls.qc_nanoCTAR_e13_5kb_all_bins.rds")


QC_Epi <- data.frame(
  sample = c(
    paste0(obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_1"]]@meta.data[["modality"]],"_nanoCTAR_e13_1"),
    paste0(obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_2"]]@meta.data[["modality"]],"_nanoCTAR_e13_2"),
    paste0(obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e15_1"]]@meta.data[["modality"]],"_nanoCTAR_e15_1"),
    paste0(obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_1"]]@meta.data[["modality"]],"_nanoCTAR_e13_1"),
    paste0(obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_2"]]@meta.data[["modality"]],"_nanoCTAR_e13_2"),
    paste0(obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e15_1"]]@meta.data[["modality"]],"_nanoCTAR_e15_1"),
    paste0(obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_1"]]@meta.data[["modality"]],"_nanoCTAR_e13_1"),
    paste0(obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_2"]]@meta.data[["modality"]],"_nanoCTAR_e13_2"),
    paste0(obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e15_1"]]@meta.data[["modality"]],"_nanoCTAR_e15_1")),
  counts = c(
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_1"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_2"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e15_1"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_1"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_2"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e15_1"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_1"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_2"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e15_1"]]@meta.data[["nCount_bins"]]
  ),
  fragment_per_feature = c(
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_1"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_2"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e15_1"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_1"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_2"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e15_1"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_1"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_2"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e15_1"]]@meta.data[["nFeature_bins"]]
  ),
  Frip = c(
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_1"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_2"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e15_1"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_1"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_2"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e15_1"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_1"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_2"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e15_1"]]@meta.data[["peak_ratio_MB"]]
  ),
  read_per_cell = c(
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_1"]]@meta.data[["total"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_2"]]@meta.data[["total"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e15_1"]]@meta.data[["total"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_1"]]@meta.data[["total"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_2"]]@meta.data[["total"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e15_1"]]@meta.data[["total"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_1"]]@meta.data[["total"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_2"]]@meta.data[["total"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e15_1"]]@meta.data[["total"]]
  ),
  duplicate = c(
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_1"]]@meta.data[["duplicate"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e13_2"]]@meta.data[["duplicate"]],
    obj.ls.qc_nanoCTAR[["ATAC_TATAGCCT_nanoCTAR_e15_1"]]@meta.data[["duplicate"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_1"]]@meta.data[["duplicate"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e13_2"]]@meta.data[["duplicate"]],
    obj.ls.qc_nanoCTAR[["H3K27ac_CCTATCCT_nanoCTAR_e15_1"]]@meta.data[["duplicate"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_1"]]@meta.data[["duplicate"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e13_2"]]@meta.data[["duplicate"]],
    obj.ls.qc_nanoCTAR[["H3K27me3_ATAGAGGC_nanoCTAR_e15_1"]]@meta.data[["duplicate"]]
  )
)

# Define the desired order
desired_order <- c("ATAC_nanoCTAR_e13_1","ATAC_nanoCTAR_e13_2","ATAC_nanoCTAR_e15_1",
                   "H3K27ac_nanoCTAR_e13_1","H3K27ac_nanoCTAR_e13_2","H3K27ac_nanoCTAR_e15_1",
                   "H3K27me3_nanoCTAR_e13_1","H3K27me3_nanoCTAR_e13_2","H3K27me3_nanoCTAR_e15_1")

# Convert the sample column to a factor with the specified order
QC_Epi$sample <- factor(QC_Epi$sample, levels = desired_order)
QC_Epi$dup_rate<-QC_Epi$duplicate/QC_Epi$read_per_cell
QC_Epi$feature_to_fragemnt<-QC_Epi$counts/QC_Epi$fragment_per_feature

p1=ggplot(QC_Epi) +
  aes(x = sample, y = counts, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)+
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 0)), 
               vjust = -5, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Number of fragments per cell") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_e13_1","ATAC_nanoCTAR_e13_2","ATAC_nanoCTAR_e15_1",
                              "H3K27ac_nanoCTAR_e13_1","H3K27ac_nanoCTAR_e13_2","H3K27ac_nanoCTAR_e15_1",
                              "H3K27me3_nanoCTAR_e13_1","H3K27me3_nanoCTAR_e13_2","H3K27me3_nanoCTAR_e15_1"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  ylim(0,5000)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p1
ggsave("/proj/user/mattia/Embryo_2/fragment_counts_plot_nanoCTAR.png", plot = p1, width = 20, height = 10, dpi = 300, units = "in")



p2=ggplot(QC_Epi) +
  aes(x = sample, y = feature_to_fragemnt, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -8, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Fragment/Features") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_e13_1","ATAC_nanoCTAR_e13_2","ATAC_nanoCTAR_e15_1",
                              "H3K27ac_nanoCTAR_e13_1","H3K27ac_nanoCTAR_e13_2","H3K27ac_nanoCTAR_e15_1",
                              "H3K27me3_nanoCTAR_e13_1","H3K27me3_nanoCTAR_e13_2","H3K27me3_nanoCTAR_e15_1"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,1)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p2
ggsave("/proj/user/mattia/Embryo_2/Fragment_to_Features_plot_nanoCTAR.png", plot = p2, width = 20, height = 10, dpi = 300, units = "in")


p3=ggplot(QC_Epi) +
  aes(x = sample, y = Frip, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)  +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -8, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Frip") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_e13_1","ATAC_nanoCTAR_e13_2","ATAC_nanoCTAR_e15_1",
                              "H3K27ac_nanoCTAR_e13_1","H3K27ac_nanoCTAR_e13_2","H3K27ac_nanoCTAR_e15_1",
                              "H3K27me3_nanoCTAR_e13_1","H3K27me3_nanoCTAR_e13_2","H3K27me3_nanoCTAR_e15_1"))+
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  ylim(0,1)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p3
ggsave("/proj/user/mattia/Embryo_2/Frip_plot_nanoCTAR.png", plot = p3, width = 20, height = 10, dpi = 300, units = "in")

p4=ggplot(QC_Epi) +
  aes(x = sample, y = read_per_cell, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)  +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -8, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("K reads per cell") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_e13_1","ATAC_nanoCTAR_e13_2","ATAC_nanoCTAR_e15_1",
                              "H3K27ac_nanoCTAR_e13_1","H3K27ac_nanoCTAR_e13_2","H3K27ac_nanoCTAR_e15_1",
                              "H3K27me3_nanoCTAR_e13_1","H3K27me3_nanoCTAR_e13_2","H3K27me3_nanoCTAR_e15_1"))+
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
ggsave("/proj/user/mattia/Embryo_2/reads_per_cell_plot_nanoCTAR.png", plot = p4, width = 20, height = 10, dpi = 300, units = "in")


p5=ggplot(QC_Epi) +
  aes(x = sample, y = dup_rate, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)  +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -3, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33","#00CC33","#00CC33")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Duplication rate per cell") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_e13_1","ATAC_nanoCTAR_e13_2","ATAC_nanoCTAR_e15_1",
                              "H3K27ac_nanoCTAR_e13_1","H3K27ac_nanoCTAR_e13_2","H3K27ac_nanoCTAR_e15_1",
                              "H3K27me3_nanoCTAR_e13_1","H3K27me3_nanoCTAR_e13_2","H3K27me3_nanoCTAR_e15_1"))+
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
ggsave("/proj/user/mattia/Embryo_2/dup_rate_plot_nanoCTAR.png", plot = p5, width = 20, height = 10, dpi = 300, units = "in")

for (experiment in names(obj.ls.qc_nanoCTAR)) {
  
  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  
  # change to UCSC style since the data was mapped to hg19
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "mm10"
  
  # add the gene information to the object
  Annotation(obj.ls.qc_nanoCTAR[[experiment]]) <- annotations
  
  obj.ls.qc_nanoCTAR[[experiment]] <- TSSEnrichment(obj.ls.qc_nanoCTAR[[experiment]], fast = FALSE)
  
  obj.ls.qc_nanoCTAR[[experiment]]  <- NucleosomeSignal(object = obj.ls.qc_nanoCTAR[[experiment]])
  
}

for (experiment in names(obj.ls.qc_nanoCTAR)) {
  obj <- obj.ls.qc_nanoCTAR[[experiment]]
  
  # Plot TSS enrichment
  p2 <- TSSPlot(obj)+ylim(0,18)

  ggsave(filename = paste0('/proj/user/mattia/Embryo_2/', experiment, '_TSS_Enrichment.png'), plot = p2)
}
