#!/usr/bin/env Rscript

# NOTE: All variables and functions are defined within the file below

source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

#***********************REMOVE MIXED POPULATION (INITIAL)**********************#

# Identify, annotate & remove mixed population clusters using UMAP of module scores
for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")){
  
  if (proj=="scRNASeq_BBN_C57B6"){
    Mixed <- list("Epithelial" = c(14,17,18,22,24,28,31,32), #28820 initial; 25214 final
                  "Fibroblasts" = c(12,15,21,29),            #36176 initial; 32550 final
                  "Myeloid" = c(20,26,27),                   #30340 initial; 29206 final
                  "Lymphoid" = c(17,21,22),                  #20661 initial; 19870 final
                  "Endothelial" = c(9,11,13,15,16,17,18))    # 6933 initial;  5786 final
  }
  
  if (proj=="scRNASeq_BBN_Rag"){
    Mixed <- list("Epithelial" = c(26,27),                                 #57184 initial; 56602 final
                  "Fibroblasts" = c(0,3,6,11,12,22),                       #17954 initial; 12618 final
                  "Myeloid" = c(0,16,19,20,25,26),                         #23171 initial; 19177 final
                  "Lymphoid" = c(),                                        #    0 initial;     0 final
                  "Endothelial" = c(0,8,12,13,14,15,18,19,20,21,22,26,27)) # 9070 initial;  5746 final
  }
  
  if (proj=="scRNASeq_Chen"){
    Mixed <- list("Epithelial" = c(9,11,26,30),              #50475 initial; 45576 final
                  "Fibroblasts" = c(13,15,18),               # 5732 initial;  5441 final
                  "Myeloid" = c(13,21,24),                   # 7475 initial;  7127 final
                  "Lymphoid" = c(22,26),                     #29420 initial; 28788 final
                  "Endothelial" = c(11,15,17,18))            #13474 initial; 12489 final
  }
  
  if (proj=="scRNASeq_GSE164557"){
    Mixed <- list("Epithelial" = c(26,29),                    #29031 initial; 28658 final
                  "Fibroblasts" = c(8,12,16,18),              # 8153 initial;  7225 final
                  "Myeloid" = c(2,3,4,6,11),                  # 1019 initial;   627 final
                  "Lymphoid" = c(4),                          #  592 initial;   536 final
                  "Endothelial" = c())                        #      initial;      final
  }
  
  if (proj=="scRNASeq_GSE132042"){
    Mixed <- list("Epithelial" = c(),                         #  4918 initial;  4918 final
                  "Fibroblasts" = c(),                        #  2885 initial;  2885 final
                  "Myeloid" = c(),                            #   184 initial;   184 final
                  "Lymphoid" = c(),                           #     0 initial;     0 final
                  "Endothelial" = c())                        #     ? initial;     ? final
    
  }
  
  if (proj=="scRNASeq_NA13_CDH12_C57B6"){
    Mixed <- list("Epithelial" = c(),                         #250009 initial; 250009 final
                  "Fibroblasts" = c(),                        #   692 initial;    692 final
                  "Myeloid" = c(),                            #  3635 initial;   3635 final
                  "Lymphoid" = c(),                           #  1315 initial;   1315 final
                  "Endothelial" = c())                        #   361 initial;    361 final
  }
  
  # Set resolution of analysis
  res <- "integrated_snn_res.1.4"
  remove_mixed(res, Mixed)
}

#******************************************************************************#
#       STEP 6: PERFORM A FINAL SUBTYPE ANALYSIS ON REQUIRED CELL TYPES        #
#******************************************************************************#

for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")){
  
  # Load the integrated seurat object
  base_data <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                              dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  res <- "integrated_snn_res.1.4"
  filt_data <- prep_data(base_data, celltype)
  sct_data <- sctransform_data(filt_data)
  plot_pre_integration(sct_data)
  integ_data <- integrate_data(sct_data, filt_data)
  integ_data <- cluster_data(integ_data)
  integ_data <- clean_data(integ_data, res)
  save_data(integ_data)
  plot_post_integration(integ_data, res)
  plot_modules(res, "All Markers")
}

#******************************************************************************#