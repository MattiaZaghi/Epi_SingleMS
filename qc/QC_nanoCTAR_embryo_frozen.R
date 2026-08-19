library(ggplot2)
library(tidyverse)
library(dplyr)
library(EnsDb.Mmusculus.v79)

obj.ls.qc_embryo<-readRDS("/date/gcb/gcb_MZ/nanoCTAR_YD/mouse/obj.ls.qc")



QC_Epi <- data.frame(
  sample = c(
    paste0(obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_sorted"]]@meta.data[["modality"]],"_embryo_sorted"),
    paste0(obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_unsorted"]]@meta.data[["modality"]],"_embryo_unsorted"),
    paste0(obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_sorted"]]@meta.data[["modality"]],"_embryo_sorted"),
    paste0(obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_unsorted"]]@meta.data[["modality"]],"_embryo_unsorted")),
  counts = c(
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_sorted"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_unsorted"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_sorted"]]@meta.data[["nCount_bins"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_unsorted"]]@meta.data[["nCount_bins"]]
    
  ),
  fragment_per_feature = c(
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_sorted"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_unsorted"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_sorted"]]@meta.data[["nFeature_bins"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_unsorted"]]@meta.data[["nFeature_bins"]]
  ),
  Frip = c(
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_sorted"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_unsorted"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_sorted"]]@meta.data[["peak_ratio_MB"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_unsorted"]]@meta.data[["peak_ratio_MB"]]
  ),
  read_per_cell = c(
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_sorted"]]@meta.data[["total"]],
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_unsorted"]]@meta.data[["total"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_sorted"]]@meta.data[["total"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_unsorted"]]@meta.data[["total"]]
  ),
  duplicate = c(
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_sorted"]]@meta.data[["duplicate"]],
    obj.ls.qc_embryo[["H3K4me3_CCTATCCT_embryo_unsorted"]]@meta.data[["duplicate"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_sorted"]]@meta.data[["duplicate"]],
    obj.ls.qc_embryo[["H3K27me3_ATAGAGGC_embryo_unsorted"]]@meta.data[["duplicate"]]
  )
)

# Define the desired order
desired_order <- c("H3K4me3_embryo_sorted","H3K4me3_embryo_unsorted",
                   "H3K27me3_embryo_sorted","H3K27me3_embryo_unsorted")

# Convert the sample column to a factor with the specified order
QC_Epi$sample <- factor(QC_Epi$sample, levels = desired_order)
QC_Epi$dup_rate<-QC_Epi$duplicate/QC_Epi$read_per_cell
QC_Epi$feature_to_fragment<-QC_Epi$counts/QC_Epi$fragment_per_feature

p1=ggplot(QC_Epi) +
  aes(x = sample, y = counts, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)+
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 0)), 
               vjust = -5, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF"))+
  scale_color_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Number of fragments per cell") +
  scale_x_discrete(labels = c("H3K4me3_embryo_sorted","H3K4me3_embryo_unsorted",
                              "H3K27me3_embryo_sorted","H3K27me3_embryo_unsorted"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,5000)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p1
ggsave("/date/gcb/gcb_MZ/nanoCTAR_YD/mouse/fragments.png", plot = p1, width = 22, height = 10, dpi = 300, units = "in")



p2=ggplot(QC_Epi) +
  aes(x = sample, y = feature_to_fragment, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)+
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -5, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF"))+
  scale_color_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Feature/Fragments") +
  scale_x_discrete(labels = c("H3K4me3_embryo_sorted","H3K4me3_embryo_unsorted",
                              "H3K27me3_embryo_sorted","H3K27me3_embryo_unsorted"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,5000)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p2
ggsave("/date/gcb/gcb_MZ/nanoCTAR_YD/mouse/feature.png", plot = p2, width = 22, height = 10, dpi = 300, units = "in")


p3=ggplot(QC_Epi) +
  aes(x = sample, y = Frip, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)+
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -5, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF"))+
  scale_color_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Frip") +
  scale_x_discrete(labels = c("H3K4me3_embryo_sorted","H3K4me3_embryo_unsorted",
                              "H3K27me3_embryo_sorted","H3K27me3_embryo_unsorted"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,5000)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p3
ggsave("/date/gcb/gcb_MZ/nanoCTAR_YD/mouse/frip.png", plot = p3, width = 22, height = 10, dpi = 300, units = "in")

p4=ggplot(QC_Epi) +
  aes(x = sample, y = read_per_cell, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)+
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 0)), 
               vjust = -5, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF"))+
  scale_color_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Reads per cell") +
  scale_x_discrete(labels = c("H3K4me3_embryo_sorted","H3K4me3_embryo_unsorted",
                              "H3K27me3_embryo_sorted","H3K27me3_embryo_unsorted"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,5000)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p4
ggsave("/date/gcb/gcb_MZ/nanoCTAR_YD/mouse/read_per_cell.png", plot = p4, width = 22, height = 10, dpi = 300, units = "in")


p5=ggplot(QC_Epi) +
  aes(x = sample, y = dup_rate, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA)+
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -5, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF"))+
  scale_color_manual(values = c("#FF0000","#FF0000", "#3399FF", "#3399FF")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Duplication rate") +
  scale_x_discrete(labels = c("H3K4me3_embryo_sorted","H3K4me3_embryo_unsorted",
                              "H3K27me3_embryo_sorted","H3K27me3_embryo_unsorted"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,5000)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p5
ggsave("/date/gcb/gcb_MZ/nanoCTAR_YD/mouse/dup_rate.png", plot = p5, width = 22, height = 10, dpi = 300, units = "in")

for (experiment in names(obj.ls.qc_nanoCTAR_K4)) {
  
  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  
  # change to UCSC style since the data was mapped to hg19
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "mm10"
  
  # add the gene information to the object
  Annotation(obj.ls.qc_nanoCTAR_K4[[experiment]]) <- annotations
  
  obj.ls.qc_nanoCTAR_K4[[experiment]] <- TSSEnrichment(obj.ls.qc_nanoCTAR_K4[[experiment]], fast = FALSE)
  
  obj.ls.qc_nanoCTAR_K4[[experiment]]  <- NucleosomeSignal(object = obj.ls.qc_nanoCTAR_K4[[experiment]])
  
}

for (experiment in names(obj.ls.qc_nanoCTAR_K4)) {
  obj <- obj.ls.qc_nanoCTAR_K4[[experiment]]
  
  # Plot TSS enrichment
  p2 <- TSSPlot(obj)+ylim(0,18)
  
  ggsave(filename = paste0('/proj/user/mattia/Embryo_2/', experiment, '_TSS_Enrichment.png'), plot = p2)
}
