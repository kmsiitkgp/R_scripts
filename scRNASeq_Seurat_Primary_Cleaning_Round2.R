#!/usr/bin/env Rscript

# NOTE: All variables and functions are defined within the file below

source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

#***********************REMOVE MIXED POPULATION (FINAL)************************#

annotations <- get_annotations(species)

# Identify, annotate & remove mixed population clusters using UMAP of module scores
for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")){
  
  if (proj=="scRNASeq_BBN_C57B6"){
    Mixed <- list("Epithelial" = c(),             #25214 initial; 25214 final
                  "Fibroblasts" = c(),            #32550 initial; 32550 final
                  "Myeloid" = c(),                #29206 initial; 29206 final
                  "Lymphoid" = c(),               #19870 initial; 19870 final
                  "Endothelial" = c())            # 5786 initial;  5786 final
  }
  if (proj=="scRNASeq_BBN_Rag"){
    Mixed <- list("Epithelial" = c(30,34),        #56602 initial; 55808 final
                  "Fibroblasts" = c(15,18),       #12618 initial; 11850 final
                  "Myeloid" = c(17),              #19177 initial; 18668 final
                  "Lymphoid" = c(),               #    0 initial;     0 final
                  "Endothelial" = c(0,14,15))     # 5746 initial; 4820 final
  }
  if (proj=="scRNASeq_Chen"){
    Mixed <- list("Epithelial" = c(8, 27,30),     #45576 initial; 43099 final
                  "Fibroblasts" = c(),            # 5441 initial;  5441 final
                  "Myeloid" = c(),                # 7127 initial;  7127 final
                  "Lymphoid" = c(),               #28788 initial; 28788 final
                  "Endothelial" = c())            #12489 initial; 12489 final
  }
  if (proj=="scRNASeq_GSE164557"){
    Mixed <- list("Epithelial" = c(),               # 28658 initial; 28658 final
                  "Fibroblasts" = c(11),            #  7225 initial;  6896 final
                  "Myeloid" = c(),                  #   627 initial;   627 final
                  "Lymphoid" = c(),                 #   536 initial;   536 final
                  "Endothelial" = c())              # initial;  final
  }
  if (proj=="scRNASeq_GSE132042"){
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
  identify_lineage(celltype)
  get_markers(annotations, res)
  plot_modules(res, celltype)
  # results_path <- seurat_results
  # prep_DESeq2(celltype)
}