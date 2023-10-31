#!/usr/bin/env Rscript

# NOTE: All variables and functions are defined within the file below

source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
#******************************************************************************#
#                          STEP 4: CLUSTER ANNOTATION                          #
#******************************************************************************#

# Load the integrated seurat object
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))

#**************************STEP 4A: CLUSTER RE-NAMING**************************#

# NOTE: If you assign same cell type to multiple clusters, verify that the
# clusters indeed are similar. This can be done by checking the clustering at a
# lower resolution. At lower resolution, if the clusters you assigned the same
# cell type are merged, then, your cell type identification is correct.
# In this example, cluster 15 and cluster 25 at resolution 1.4 seemed similar.
# Indeed, the 2 clusters are merged into cluster 5 at resolution 1.2 or
# cluster 0 at resolution 1

# (i) DO NOT USE symbols like "/" in cluster names as it will create errors when
# saving files based on cluster names.
# (ii) You can use Seurat::RenameIdents() and then add cell type identity to
# metadata using integrated_seurat$cell_type <- Idents(integrated_seurat1).

# Easier alternative is to use custom code below which has multiple benefits:
# (a) Using RenameIdents() permanently changes Idents. So, if we want to change
# annotation, we have to read the Seurat object again. To avoid this, we can
# store the annotated data in new object "integrated_seurat1" but this increases
# memory used by R
# (b) If we save the annotated "integrated_seurat1" as rds, it will take another
# 10GB space in cluster.
# By using custom code, we aren't altering the idents of the seurat object.
# This allows to rename the clusters freely on the same seurat object.
# Also, we can save simply overwrite "integrated_seurat" as rds saving space

# Renaming using RenameIdents() is NOT recommended
# integrated_seurat <- SeuratObject::RenameIdents(object = integrated_seurat,
#                                    "0" = "")

if (proj=="scRNASeq_BBN_C57B6"){
  clusters <- list("Fibroblasts" = c(3,5,6,13,16,17,23,26,29,30,33,35,42,54),
                   "Epithelial" = c(0,8,11,14,15,18,25,34,39,51),
                   "Myeloid - Macrophages, DCs" = c(9,10,20,32,37,44,45,49),
                   "Myeloid - MDSC" = c(4,7,12,21),
                   "Myeloid - Mast" = c(),
                   "Lymphoid - B" = c(2,47,48),
                   "Lymphoid - T" = c(19,22,24,27,28,41),
                   "Lymphoid - NK" = c(31),
                   "Endothelial" = c(1,46,50),
                   "Endothelial - Lymphatic" = c(38),
                   "Neurons" = c(43),
                   "Muscle" = c(36),
                   "Unclassified" = c(40,52,53))
}

if (proj=="scRNASeq_BBN_Rag"){
  clusters <- list("Fibroblasts" = c(7,8,12,22,26,31,40,43),
                   "Epithelial" = c(2,3,4,5,6,9,10,11,13,15,16,17,19,21,27,29,32,36,42,37),
                   "Myeloid - Macrophages, DCs" = c(0,18,28,33,34,41,20,45),
                   "Myeloid - MDSC" = c(1,23),
                   "Myeloid - Mast" = c(),
                   "Lymphoid - B" = c(),
                   "Lymphoid - T" = c(),
                   "Lymphoid - NK" = c(),
                   "Endothelial" = c(14,30,35,39,24,38),
                   "Endothelial - Lymphatic" = c(44),
                   "Neurons" = c(),
                   "Muscle" = c(25),
                   "Unclassified" = c())
}

if (proj=="scRNASeq_GSE164557"){
  clusters <- list("Fibroblasts" = c(3,4,12,14,29),
                   "Epithelial" = c(0,1,2,5,6,7,8,9,10,11,13,15,16,17,18,19,24,25,28,30),
                   "Myeloid - Macrophages, DCs" = c(21,22,23),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Mast" = c(),
                   "Lymphoid - B" = c(),
                   "Lymphoid - T" = c(20),
                   "Lymphoid - NK" = c(),
                   "Endothelial" = c(26),
                   "Endothelial - Lymphatic" = c(),
                   "Neurons" = c(),
                   "Muscle" = c(27),
                   "Unclassified" = c())
}

if (proj=="scRNASeq_NA13_CDH12_C57B6"){
  clusters <- list("Fibroblasts" = c(32),
                   "Epithelial" = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
                                    18,19,20,21,22,23,24,25,27,30,31,33,34,36,37,
                                    39,78),
                   "Myeloid - Macrophages, DCs" = c(26, 28),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Mast" = c(),
                   "Lymphoid - B" = c(),
                   "Lymphoid - T" = c(29),
                   "Lymphoid - NK" = c(),
                   "Endothelial" = c(35),
                   "Endothelial - Lymphatic" = c(),
                   "Neurons" = c(),
                   "Muscle" = c(38),
                   "Unclassified" = c())
}

if (proj=="scRNASeq_Chen"){
  clusters <- list("Fibroblasts" = c(22,27,37,38),
                   "Epithelial" = c(0,2,6,7,8,12,14,15,18,23,24,28,30,40,46,47,5),
                   "Myeloid - Macrophages, DCs" = c(19,29,33),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Mast" = c(20),
                   "Lymphoid - B" = c(16,36),
                   "Lymphoid - T" = c(3,4,9,10,13,31,35,41),
                   "Lymphoid - NK" = c(25),
                   "Endothelial" = c(1,11,26,34,39),
                   "Endothelial - Lymphatic" = c(43),
                   "Neurons" = c(44),
                   "Muscle" = c(17,21,32,42),
                   "Unclassified" = c(45))
}

# if (proj=="scRNASeq_GSE132042"){
#   clusters <- list("Fibroblasts" = c(5,6,7,8,10,13),
#                    "Epithelial" = c(0,1,2,3,4,9),
#                    "Myeloid - Macrophages, DCs" = c(11,16),
#                    "Myeloid - MDSC" = c(),
#                    "Myeloid - Mast" = c(),
#                    "Lymphoid - B" = c(),
#                    "Lymphoid - T" = c(),
#                    "Lymphoid - NK" = c(),
#                    "Endothelial" = c(14,15),
#                    "Endothelial - Lymphatic" = c(),
#                    "Neurons" = c(),
#                    "Muscle" = c(12),
#                    "Unclassified" = c())
# }

# Make sure you have assigned all clusters to one of the cell types
# NOTE: "integrated_snn_res.1.4" etc are stored as factors. 
# So, use as.character() and then as.numeric() to get accurate cluster values
list_1 <- integrated_seurat@meta.data %>% 
  dplyr::count(integrated_snn_res.1.4) %>% 
  dplyr::select(identity(1)) %>% 
  unlist(use.names=FALSE) %>% 
  as.character() %>% 
  as.numeric() %>% 
  sort()

list_2 <- clusters %>% 
  unlist(., use.names=FALSE) %>% 
  sort()

# Proceed with annotation ONLY if all clusters have been renamed
if (identical(list_1, list_2)){
  print("All Clusters have been annotated")
  
  # Extract metadata from Seurat object, assign appropriate resolution to
  # seurat_clusters column and add cell_type, sub_type columns
  data <- integrated_seurat@meta.data %>% 
    dplyr::mutate(seurat_clusters = !!rlang::sym("integrated_snn_res.1.4"),
                  cell_type = NA, sub_type = NA, cell_class = NA)
  
  # Assign cell type based on cluster numbers within seurat_clusters column
  for (j in 1:nrow(data)){
    for (i in 1:length(clusters)){
      if (as.numeric(as.character(data$seurat_clusters[j])) %in% clusters[[i]]){
        data$cell_type[j] <- names(clusters[i])
      }
    }
  }
  
  # Check summary of cell counts
  print(data %>% dplyr::count(cell_type) %>% dplyr::arrange(n))
  cat("\n")
  
  # Import the metadata into Seurat object and save it
  integrated_seurat@meta.data <- data
  saveRDS(integrated_seurat, paste0(seurat_results, "integrated_seurat_snn.rds"))
} else {
  cat("\nYou missed annotating these clusters:\t", setdiff(list_1, list_2))
  cat("\nThese clusters are not present in data:\t", setdiff(list_2, list_1))
}

#******************************************************************************#
#       STEP 5: PERFORM INITIAL SUBTYPE ANALYSIS ON REQUIRED CELL TYPES        #
#******************************************************************************#

# NOTE: We have defined "Myeloid - MDSCs, Myeloid - Mast Cells etc..so, use
# "Myeloid" in for loop, NOT "Myeloid Cells"

# Load the integrated seurat object
base_data <- readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))

for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")){
  
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