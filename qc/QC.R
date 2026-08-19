library(ggplot2)
library(tidyverse)
library(dplyr)

#load all the datasets to create one unique table to compare different QC Metric from the different datasets

combined.obj.ls_H3K27ac_all_mod_PT<-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_H3K27ac_all_mod_paired_tag.rds")

combined.obj.ls_H3K27me3_all_mod_PT<-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_H3K27me3_all_mod_paired_tag.rds")

combined.obj.ls_all_mod_multiome<-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_all_mod_embryo_10x.rds")

combined.obj.ls_all_mod_nanoCT<-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_all_mod_nanoCT.rds")

combined.obj.ls_all_mod_embryo_old <-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_all_amb.rds")

combined.obj.ls_scCUTTAG_Pro_H3K27ac<-readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_epi_scCut&Tag_pro_H3K27ac_5kb_bins.rds")

combined.obj.ls_scCUTTAG_Pro <- readRDS("/date/gcb/gcb_MZ/multiNanoCT/samples/Analysis/combined.obj.ls_epi_scCut&Tag_pro_H3K4me3_27me3_5kb_bins.rds")

combined.obj.ls_scCUTTAG_Pro[["H3K27ac_CCTATCCT"]]<-combined.obj.ls_scCUTTAG_Pro_H3K27ac[["H3K27ac_CCTATCCT"]]

combined.obj.ls_new_nanoCTAR<-readRDS("/proj/user/mattia/Embryo_2/combined.obj.ls_epi_mod_nanoCTAR_e13_5kb_bins_epi_all.rds")
combined.obj.ls_new_nanoCTAR_K4<-readRDS("/proj/user/mattia/Embryo_2/combined.obj.ls_epi_mod_nanoCTAR_e13_K4_5kb_bins.rds")
combined.obj.ls_new_nanoCTR_K4<-readRDS("/proj/user/mattia/Embryo_2/combined.obj.ls_epi_mod_nanoCTR_e13_K4_5kb_bins.rds")



QC_Epi <- data.frame(
  sample = c(
    paste0(combined.obj.ls_new_nanoCTAR[["ATAC_TATAGCCT"]]@meta.data[["modality"]], "_nanoCTAR_new"),
    paste0(combined.obj.ls_new_nanoCTAR_K4[["ATAC_TATAGCCT"]]@meta.data[["modality"]], "_nanoCTAR_new"),
    paste0(combined.obj.ls_all_mod_embryo_old[["ATAC_TATAGCCT"]]@meta.data[["modality"]], "_nanoCTAR_old"),
    paste0(combined.obj.ls_all_mod_nanoCT[["ATAC_TATAGCCT"]]@meta.data[["modality"]], "_nanoCT"),
    paste0(combined.obj.ls_all_mod_multiome[["ATAC_TATAGCCT"]]@meta.data[["modality"]], "_Embryo_10x_multiome"),
    paste0(combined.obj.ls_new_nanoCTAR[["H3K27ac_CCTATCCT"]]@meta.data[["modality"]], "_nanoCTAR_new"),
    paste0(combined.obj.ls_all_mod_embryo_old[["H3K27ac_CCTATCCT"]]@meta.data[["modality"]], "_nanoCTAR_old"),
    paste0(combined.obj.ls_all_mod_nanoCT[["H3K27ac_ATAGAGGC"]]@meta.data[["modality"]], "_nanoCT"),
    paste0(combined.obj.ls_H3K27ac_all_mod_PT[["H3K27ac_CCTATCCT"]]@meta.data[["modality"]], "_Droplet_Pair-Tag"),
    paste0(combined.obj.ls_scCUTTAG_Pro[["H3K27ac_CCTATCCT"]]@meta.data[["modality"]], "_scCut&Tag-pro"),
    paste0(combined.obj.ls_new_nanoCTAR[["H3K27me3_ATAGAGGC"]]@meta.data[["modality"]], "_nanoCTAR_new"),
    paste0(combined.obj.ls_new_nanoCTAR_K4[["H3K27me3_ATAGAGGC"]]@meta.data[["modality"]], "_nanoCTAR_new"),
    paste0(combined.obj.ls_all_mod_embryo_old[["H3K27me3_ATAGAGGC"]]@meta.data[["modality"]], "_nanoCTAR_old"),
    paste0(combined.obj.ls_all_mod_nanoCT[["H3K27me3_CCTATCCT"]]@meta.data[["modality"]], "_nanoCT"),
    paste0(combined.obj.ls_H3K27me3_all_mod_PT[["H3K27me3_ATAGAGGC"]]@meta.data[["modality"]], "_Droplet_Pair-Tag"),
    paste0(combined.obj.ls_scCUTTAG_Pro[["H3K27me3_TATAGCCT_PBMCs"]]@meta.data[["modality"]], "_scCut&Tag-pro"),
    paste0(combined.obj.ls_new_nanoCTAR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["modality"]], "_nanoCTAR_new"),
    paste0(combined.obj.ls_new_nanoCTR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["modality"]], "_nanoCTAR_new"),
    paste0(combined.obj.ls_scCUTTAG_Pro[["H3K4me3_ATAGAGGC_PBMCs"]]@meta.data[["modality"]], "_scCut&Tag-pro")
  ),
  counts = c(
    combined.obj.ls_new_nanoCTAR[["ATAC_TATAGCCT"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_new_nanoCTAR_K4[["ATAC_TATAGCCT"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_all_mod_embryo_old[["ATAC_TATAGCCT"]]@meta.data[["nCount_ATAC_bins"]],
    combined.obj.ls_all_mod_nanoCT[["ATAC_TATAGCCT"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_all_mod_multiome[["ATAC_TATAGCCT"]]@meta.data[["nCount_ATAC_bins"]],
    combined.obj.ls_new_nanoCTAR[["H3K27ac_CCTATCCT"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27ac_CCTATCCT"]]@meta.data[["nCount_H3K27ac_bins"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27ac_ATAGAGGC"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_H3K27ac_all_mod_PT[["H3K27ac_CCTATCCT"]]@meta.data[["nCount_H3K27ac_bins"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27ac_CCTATCCT"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_new_nanoCTAR[["H3K27me3_ATAGAGGC"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K27me3_ATAGAGGC"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27me3_ATAGAGGC"]]@meta.data[["nCount_H3K27me3_bins"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27me3_CCTATCCT"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_H3K27me3_all_mod_PT[["H3K27me3_ATAGAGGC"]]@meta.data[["nCount_H3K27me3_bins"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27me3_TATAGCCT_PBMCs"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_new_nanoCTR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["nCount_bins"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K4me3_ATAGAGGC_PBMCs"]]@meta.data[["nCount_bins"]]
  ),
  fragment_per_feature = c(
    combined.obj.ls_new_nanoCTAR[["ATAC_TATAGCCT"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_new_nanoCTAR_K4[["ATAC_TATAGCCT"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_all_mod_embryo_old[["ATAC_TATAGCCT"]]@meta.data[["nFeature_ATAC_bins"]],
    combined.obj.ls_all_mod_nanoCT[["ATAC_TATAGCCT"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_all_mod_multiome[["ATAC_TATAGCCT"]]@meta.data[["nFeature_ATAC_bins"]],
    combined.obj.ls_new_nanoCTAR[["H3K27ac_CCTATCCT"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27ac_CCTATCCT"]]@meta.data[["nFeature_H3K27ac_bins"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27ac_ATAGAGGC"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_H3K27ac_all_mod_PT[["H3K27ac_CCTATCCT"]]@meta.data[["nFeature_H3K27ac_bins"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27ac_CCTATCCT"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_new_nanoCTAR[["H3K27me3_ATAGAGGC"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K27me3_ATAGAGGC"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27me3_ATAGAGGC"]]@meta.data[["nFeature_H3K27me3_bins"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27me3_CCTATCCT"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_H3K27me3_all_mod_PT[["H3K27me3_ATAGAGGC"]]@meta.data[["nFeature_H3K27me3_bins"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27me3_TATAGCCT_PBMCs"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_new_nanoCTR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["nFeature_bins"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K4me3_ATAGAGGC_PBMCs"]]@meta.data[["nFeature_bins"]]
  ),
  Frip = c(
    combined.obj.ls_new_nanoCTAR[["ATAC_TATAGCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_new_nanoCTAR_K4[["ATAC_TATAGCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_all_mod_embryo_old[["ATAC_TATAGCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_all_mod_nanoCT[["ATAC_TATAGCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_all_mod_multiome[["ATAC_TATAGCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_new_nanoCTAR[["H3K27ac_CCTATCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27ac_CCTATCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27ac_ATAGAGGC"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_H3K27ac_all_mod_PT[["H3K27ac_CCTATCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27ac_CCTATCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_new_nanoCTAR[["H3K27me3_ATAGAGGC"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K27me3_ATAGAGGC"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27me3_ATAGAGGC"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27me3_CCTATCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_H3K27me3_all_mod_PT[["H3K27me3_ATAGAGGC"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27me3_TATAGCCT_PBMCs"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_new_nanoCTR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["peak_ratio_MB"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K4me3_ATAGAGGC_PBMCs"]]@meta.data[["peak_ratio_MB"]]
  ),
  read_per_cell = c(
    combined.obj.ls_new_nanoCTAR[["ATAC_TATAGCCT"]]@meta.data[["total"]],
    combined.obj.ls_new_nanoCTAR_K4[["ATAC_TATAGCCT"]]@meta.data[["total"]],
    combined.obj.ls_all_mod_embryo_old[["ATAC_TATAGCCT"]]@meta.data[["total"]],
    combined.obj.ls_all_mod_nanoCT[["ATAC_TATAGCCT"]]@meta.data[["total"]],
    combined.obj.ls_all_mod_multiome[["ATAC_TATAGCCT"]]@meta.data[["total"]],
    combined.obj.ls_new_nanoCTAR[["H3K27ac_CCTATCCT"]]@meta.data[["total"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27ac_CCTATCCT"]]@meta.data[["total"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27ac_ATAGAGGC"]]@meta.data[["total"]],
    combined.obj.ls_H3K27ac_all_mod_PT[["H3K27ac_CCTATCCT"]]@meta.data[["total"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27ac_CCTATCCT"]]@meta.data[["total"]],
    combined.obj.ls_new_nanoCTAR[["H3K27me3_ATAGAGGC"]]@meta.data[["total"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K27me3_ATAGAGGC"]]@meta.data[["total"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27me3_ATAGAGGC"]]@meta.data[["total"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27me3_CCTATCCT"]]@meta.data[["total"]],
    combined.obj.ls_H3K27me3_all_mod_PT[["H3K27me3_ATAGAGGC"]]@meta.data[["total"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27me3_TATAGCCT_PBMCs"]]@meta.data[["total"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["total"]],
    combined.obj.ls_new_nanoCTR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["total"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K4me3_ATAGAGGC_PBMCs"]]@meta.data[["total"]]
  ),
  duplicate = c(
    combined.obj.ls_new_nanoCTAR[["ATAC_TATAGCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_new_nanoCTAR_K4[["ATAC_TATAGCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_all_mod_embryo_old[["ATAC_TATAGCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_all_mod_nanoCT[["ATAC_TATAGCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_all_mod_multiome[["ATAC_TATAGCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_new_nanoCTAR[["H3K27ac_CCTATCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27ac_CCTATCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27ac_ATAGAGGC"]]@meta.data[["duplicate"]],
    combined.obj.ls_H3K27ac_all_mod_PT[["H3K27ac_CCTATCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27ac_CCTATCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_new_nanoCTAR[["H3K27me3_ATAGAGGC"]]@meta.data[["duplicate"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K27me3_ATAGAGGC"]]@meta.data[["duplicate"]],
    combined.obj.ls_all_mod_embryo_old[["H3K27me3_ATAGAGGC"]]@meta.data[["duplicate"]],
    combined.obj.ls_all_mod_nanoCT[["H3K27me3_CCTATCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_H3K27me3_all_mod_PT[["H3K27me3_ATAGAGGC"]]@meta.data[["duplicate"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K27me3_TATAGCCT_PBMCs"]]@meta.data[["duplicate"]],
    combined.obj.ls_new_nanoCTAR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_new_nanoCTR_K4[["H3K4me3_CCTATCCT"]]@meta.data[["duplicate"]],
    combined.obj.ls_scCUTTAG_Pro[["H3K4me3_ATAGAGGC_PBMCs"]]@meta.data[["duplicate"]]
  )
)

# Define the desired order
desired_order <- c("ATAC_nanoCTAR_new","ATAC_nanoCTAR_old","ATAC_nanoCT","ATAC_Embryo_10x_multiome","H3K27ac_nanoCTAR_new","H3K27ac_nanoCTAR_old",    
                   "H3K27ac_nanoCT","H3K27ac_Droplet_Pair-Tag","H3K27ac_scCut&Tag-pro",
                   "H3K27me3_nanoCTAR_new","H3K27me3_nanoCTAR_old","H3K27me3_nanoCT","H3K27me3_Droplet_Pair-Tag","H3K27me3_scCut&Tag-pro",
                   "H3K4me3_nanoCTAR_new","H3K4me3_scCut&Tag-pro")

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
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                               "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                                "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Number of fragments per cell") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_new","ATAC_nanoCTAR_old","ATAC_nanoCT","ATAC_Embryo_10x_multiome","H3K27ac_nanoCTAR_new","H3K27ac_nanoCTAR_old",    
                              "H3K27ac_nanoCT","H3K27ac_Droplet_Pair-Tag","H3K27ac_scCut&Tag-pro",
                              "H3K27me3_nanoCTAR_new","H3K27me3_nanoCTAR_old","H3K27me3_nanoCT","H3K27me3_Droplet_Pair-Tag","H3K27me3_scCut&Tag-pro",
                              "H3K4me3_nanoCTAR_new","H3K4me3_scCut&Tag-pro"))+
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
ggsave("/proj/user/mattia/Embryo_2/fragment_counts_plot.png", plot = p1, width = 20, height = 10, dpi = 300, units = "in")

p2=ggplot(QC_Epi) +
  aes(x = sample, y = feature_to_fragemnt, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -8, size = 5, color = "black") +
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                               "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                                "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Fragment/Features") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_new","ATAC_nanoCTAR_old","ATAC_nanoCT","ATAC_Embryo_10x_multiome","H3K27ac_nanoCTAR_new","H3K27ac_nanoCTAR_old",    
                              "H3K27ac_nanoCT","H3K27ac_Droplet_Pair-Tag","H3K27ac_scCut&Tag-pro",
                              "H3K27me3_nanoCTAR_new","H3K27me3_nanoCTAR_old","H3K27me3_nanoCT","H3K27me3_Droplet_Pair-Tag","H3K27me3_scCut&Tag-pro",
                              "H3K4me3_nanoCTAR_new","H3K4me3_scCut&Tag-pro"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  #ylim(0,1)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p2
ggsave("/proj/user/mattia/Embryo_2/Fragment_to_Features_plot.png", plot = p2, width = 20, height = 10, dpi = 300, units = "in")

p3=ggplot(QC_Epi) +
  aes(x = sample, y = Frip, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -4, size = 5, color = "black")+
  scale_fill_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                               "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  scale_color_manual(values = c("#FF0000","#FF0000", "#FF0000","#FF0000","#3399FF", "#3399FF", "#3399FF","#3399FF","#3399FF",
                                "#00CC33","#00CC33","#00CC33","#00CC33","#00CC33","#FF9933","#FF9933")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Fraction of reads into peaks") +
  scale_x_discrete(labels = c("ATAC_nanoCTAR_new","ATAC_nanoCTAR_old","ATAC_nanoCT","ATAC_Embryo_10x_multiome","H3K27ac_nanoCTAR_new","H3K27ac_nanoCTAR_old",    
                              "H3K27ac_nanoCT","H3K27ac_Droplet_Pair-Tag","H3K27ac_scCut&Tag-pro",
                              "H3K27me3_nanoCTAR_new","H3K27me3_nanoCTAR_old","H3K27me3_nanoCT","H3K27me3_Droplet_Pair-Tag","H3K27me3_scCut&Tag-pro",
                              "H3K4me3_nanoCTAR_new","H3K4me3_scCut&Tag-pro"))+
  #scale_y_continuous(breaks = seq(0, 2000, by = 1000)) +
  ylim(0,1)+
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
p3
ggsave("/proj/user/mattia/Embryo_2/Frip_plot.png", plot = p3, width = 20, height = 10, dpi = 300, units = "in")


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

for (experiment in names(obj.ls.qc_nanoCTAR_K4)) {
  
  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  
  # change to UCSC style since the data was mapped to hg19
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "mm10"
  
  # add the gene information to the object
  Annotation(obj.ls.qc_nanoCTAR_K4[[experiment]]) <- annotations
  
  obj.ls.qc_nanoCTAR_K4[[experiment]] <- TSSEnrichment(obj.ls.qc_nanoCTAR_K4[[experiment]], fast = FALSE)
}

for (experiment in names(obj.ls.qc_nanoCTAR_K4)) {
  obj <- obj.ls.qc_nanoCTAR_K4[[experiment]]
  
  # Plot TSS enrichment
  p2 <- TSSPlot(obj)
  
  ggsave(filename = paste0('/proj/user/mattia/Embryo_2/', experiment, '_TSS_Enrichment.png'), plot = p2)
}
