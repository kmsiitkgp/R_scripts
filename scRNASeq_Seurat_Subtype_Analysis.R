#!/usr/bin/env Rscript

# NOTE: All variables and functions are defined within the file below
# https://github.com/satijalab/seurat/issues/2101

source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

#******************************************************************************#
#                         STEP 10: SUBTYPE ANNOTATION                          #
#******************************************************************************#

# Use UMAP@res1.4, UMAP split by Condition, FindMarkers(), vst plots to decide subtypes
for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")){
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  if (proj == "scRNASeq_BBN_C57B6"){
    if (celltype == "Fibroblasts"){
      clusters <- list("Fibroblasts - Inflammatory" = c(13,15,16,24,26),
                       "Fibroblasts - Myofibrotic" = c(4,27),
                       "Fibroblasts - Proliferating" = c(20),
                       "Fibroblasts - CAF-A" = c(11,18,19,23,25),
                       "Fibroblasts - CAF-B" = c(5,6),
                       "Fibroblasts - nonCAF" = c(0,2,22,1,7,10,12,3,8,9,14,28,29),# based on UMAP split by "Condition", cell% in excel file
                       "Fibroblasts - Mesothelial-like" = c(17,21),
                       "Unclassified" = c())
    }
    if (celltype == "Epithelial"){
      clusters <- list("Epithelial - Basal" = c(8,10,19,20,21,22,23,4,6,7,11,12,2,17,14),
                       "Epithelial - Luminal" = c(1,3,9,13),
                       "Epithelial - Umbrella" = c(16,24),
                       "Epithelial - Inflammatory" = c(0),
                       "Epithelial - Proliferating" = c(5),
                       "Epithelial - EMT-like" = c(18),
                       "Unclassified" = c(15,25,26,27,28,29,30))
    }
    if (celltype == "Myeloid"){
      clusters <- list("Myeloid - Macrophage" = c(3,9,14,15,17,19,22,25,4,11,23),
                       "Myeloid - cDCs" = c(13,18,27,28),
                       "Myeloid - pDCs" = c(29),
                       "Myeloid - MDSC" = c(2,5,6,10,16,7,8,0,1,12,26,30),
                       "Myeloid - Mast" = c(),
                       "Unclassified" = c(20,21,24))
    }
    if (celltype == "Lymphoid"){
      clusters <- list("Lymphoid - CD4,CD8 Naive" = c(7,9),
                       "Lymphoid - CD4 Naive" = c(),
                       "Lymphoid - CD4 Helper" = c(3,15),
                       "Lymphoid - CD4 Treg" = c(1),
                       "Lymphoid - CD8 Naive" = c(),
                       "Lymphoid - CD8 Cytotoxic" = c(0,12,20),
                       "Lymphoid - Gamma Delta" = c(13),
                       "Lymphoid - NKT" = c(4,14,22),
                       "Lymphoid - NK" = c(10,11),
                       "Lymphoid - B" = c(2,5,6,8,18,21),
                       "Lymphoid - Plasma" = c(19),
                       "Unclassified" = c(16,17,23,24))
    }
    if (celltype == "Endothelial"){
      clusters <- list("Endothelial" = c(0,1,3,4,5,6,7,8,9,10,11,12),
                       "Endothelial - Lymphatic" = c(2,13,14))
    }
  }
  
  if (proj == "scRNASeq_BBN_Rag"){
    if (celltype == "Fibroblasts"){
      clusters <- list("Fibroblasts - Inflammatory" = c(13,1,5,4,11),
                       "Fibroblasts - Myofibrotic" = c(2,7),
                       "Fibroblasts - Proliferating" = c(17),
                       "Fibroblasts - CAF-A" = c(21),
                       "Fibroblasts - CAF-B" = c(0,3,9,12),
                       "Fibroblasts - CAF-C" = c(6,16,19),
                       "Fibroblasts - nonCAF" = c(), # based on UMAP split by "Condition", cell% in excel file
                       "Fibroblasts - Mesothelial-like" = c(8,10,20),
                       "Unclassified" = c(14))
    }
    if (celltype == "Epithelial"){
      clusters <- list("Epithelial - Basal" = c(0,2,3,5,7,8,9,11,16,18,20,22,24,25,26,27,14,29),
                       "Epithelial - Luminal" = c(13,21,32,6,12,19,17,1),
                       "Epithelial - Umbrella" = c(23),
                       "Epithelial - Inflammatory" = c(10),
                       "Epithelial - Proliferating" = c(4),
                       "Epithelial - EMT-like" = c(31),
                       "Unclassified" = c(15,28,33))
    }
    if (celltype == "Myeloid"){
      clusters <- list("Myeloid - Macrophage" = c(5,9,10,13,14,16,20,0,1,6,19),
                       "Myeloid - cDCs" = c(15,22),
                       "Myeloid - pDCs" = c(25),
                       "Myeloid - MDSC" = c(2,3,4,7,8,11,12,18),
                       "Myeloid - Mast" = c(),
                       "Unclassified" = c(21,23,24))
    }
    # There are no lymphoid cells in Rag mice
    if (celltype == "Endothelial"){
      clusters <- list("Endothelial" = c(1,2,3,4,5,7,8,9,10,11,12,13,16),
                       "Endothelial - Lymphatic" = c(6),
                       "Unclassified" = c(17))
    }
  }
  
  if (proj == "scRNASeq_Chen"){
    if (celltype == "Fibroblasts"){
      clusters <- list("Fibroblasts - Inflammatory" = c(),
                       "Fibroblasts - Myofibrotic" = c(7),
                       "Fibroblasts - Proliferating" = c(),
                       "Fibroblasts - CAFs" = c(0,4,2,1,6,10,8,5,11,3,9,13,16),
                       "Fibroblasts - CAF-A" = c(),
                       "Fibroblasts - CAF-B" = c(),
                       "Fibroblasts - nonCAF" = c(12,14),# based on UMAP split by "Condition", cell% in excel file
                       "Fibroblasts - Mesothelial-like" = c(),
                       "Unclassified" = c(15))
    }
    if (celltype == "Epithelial"){
      clusters <- list("Epithelial - Basal" = c(0),
                       "Epithelial - Luminal" = c(1,2,5,7,9,10,11,13,21,22,
                                                  3,14,16,23,18,20,17,
                                                  12,15,6,19,24,26),
                       "Epithelial - Umbrella" = c(),
                       "Epithelial - Inflammatory" = c(),
                       "Epithelial - Proliferating" = c(4,25),
                       "Epithelial - EMT-like" = c(),
                       "Unclassified" = c(28,29,31))
    }
    if (celltype == "Myeloid"){
      clusters <- list("Myeloid - Macrophage" = c(0,1,5,22,6,12,13,16),
                       "Myeloid - cDCs" = c(15,17),
                       "Myeloid - pDCs" = c(18),
                       "Myeloid - MDSC" = c(3,10,11,19,21),
                       "Myeloid - Mast" = c(2,4,7,8,9,14),
                       "Unclassified" = c(20,23))
    }
    if (celltype == "Lymphoid"){
      clusters <- list("Lymphoid - CD4,CD8 Naive" = c(6,21),
                       "Lymphoid - CD4 Naive" = c(),
                       "Lymphoid - CD4 Helper" = c(),
                       "Lymphoid - CD4 Treg" = c(2,8,16),
                       "Lymphoid - CD8 Naive" = c(),
                       "Lymphoid - CD8 Cytotoxic" = c(0,1,5,11,15,23),
                       "Lymphoid - Gamma Delta" = c(),
                       "Lymphoid - NKT" = c(),
                       "Lymphoid - NK" = c(17,18,24),
                       "Lymphoid - B" = c(7,12),
                       "Lymphoid - Plasma" = c(19),
                       "Lymphoid - T" = c(3,4,9,10,13,14,20,22),
                       "Unclassified" = c(25,26,27,28))
    }
    if (celltype == "Endothelial"){
      clusters <- list("Endothelial" = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17),
                       "Endothelial - Lymphatic" = c(16))
    }
  }
  
  if (proj == "scRNASeq_GSE164557"){
    if (celltype == "Fibroblasts"){
      clusters <- list("Fibroblasts - Inflammatory" = c(),
                       "Fibroblasts - Myofibrotic" = c(6,13,14),
                       "Fibroblasts - Proliferating" = c(),
                       "Fibroblasts - CAF-A" = c(),
                       "Fibroblasts - CAF-B" = c(),
                       "Fibroblasts - CAF-C" = c(),
                       "Fibroblasts - nonCAF" = c(0,1,2,3,4,5,7,8,9,10,12,15,16,17), 
                       "Fibroblasts - Mesothelial-like" = c(),
                       "Unclassified" = c())
    }
    if (celltype == "Epithelial"){
      clusters <- list("Epithelial - Basal" = c(1,2,3,5,6,7,11,13,16,19,28,29),
                       "Epithelial - Luminal" = c(0,4,8,9,10,12,14,15,17,18,20,22,23,24,25,26,31),
                       "Epithelial - Umbrella" = c(21,27),
                       "Epithelial - Inflammatory" = c(),
                       "Epithelial - Proliferating" = c(),
                       "Epithelial - EMT-like" = c(),
                       "Unclassified" = c(30))
    }
    if (celltype == "Myeloid"){
      clusters <- list("Myeloid - Macrophage" = c(0,2,4),
                       "Myeloid - cDCs" = c(1,3,5,6,7),
                       "Myeloid - pDCs" = c(),
                       "Myeloid - MDSC" = c(),
                       "Myeloid - Mast" = c(),
                       "Unclassified" = c())
    }
    if (celltype == "Lymphoid"){
      clusters <- list("Lymphoid - CD4,CD8 Naive" = c(),
                       "Lymphoid - CD4 Naive" = c(),
                       "Lymphoid - CD4 Helper" = c(3),
                       "Lymphoid - CD4 Treg" = c(4),
                       "Lymphoid - CD8 Naive" = c(),
                       "Lymphoid - CD8 Cytotoxic" = c(1),
                       "Lymphoid - Gamma Delta" = c(2),
                       "Lymphoid - NKT" = c(),
                       "Lymphoid - NK" = c(0),
                       "Lymphoid - B" = c(),
                       "Lymphoid - Plasma" = c(),
                       "Unclassified" = c())
    }
    if (celltype == "Endothelial"){
      clusters <- list("Endothelial" = c(),
                       "Endothelial - Lymphatic" = c(),
                       "Unclassified" = c())
    }
  }
  
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
      dplyr::mutate(seurat_clusters = !!rlang::sym("integrated_snn_res.1.4"))
    
    # Assign cell type based on cluster numbers within seurat_clusters column
    for (j in 1:nrow(data)){
      for (i in 1:length(clusters)){
        if (as.numeric(as.character(data$seurat_clusters[j])) %in% clusters[[i]]){
          data$sub_type[j] <- names(clusters[i])
        }
      }
    }
    
    # Identify cross labelled cells in cell_class column of metadata
    # Some cells may have cell_type "Myeloid - MDSC" but sub_type "Myeloid- cDcs"
    integrated_seurat@meta.data <- data %>%
      dplyr::mutate(cell_class = dplyr::case_when(sub_type == "Unclassified" ~ "Unclassified",
                                                  cell_type == sub_type |
                                                    cell_type == "Lymphoid - B" & grepl(pattern = "B|Plasma", x=sub_type) |
                                                    cell_type == "Lymphoid - T" & grepl(pattern = "CD4|CD8|Gamma|NKT", x=sub_type) |
                                                    cell_type == "Myeloid - Macrophages, DCs" & grepl(pattern = "Macrophage|cDCs|pDCs|MDSC", x=sub_type) |
                                                    cell_type == "Epithelial" & grepl(pattern = "Epithelial", x=sub_type) |
                                                    cell_type == "Fibroblasts" & grepl(pattern = "Fibroblasts", x=sub_type) ~ gsub(pattern = "\ -.*",replacement = "",x = cell_type),
                                                  TRUE ~ "Mixed"))
    
    # Check what has been re-annotated. Notice how some cells were wrongly
    # annotated prior to subtype  analysis as their sub_type doesn't match their
    # cell_type. This is visible clearly in myeloid & lymphoid cells.
    print(integrated_seurat@meta.data %>% dplyr::count(cell_class, cell_type, sub_type))
    
    # Import the metadata into Seurat object and save it
    saveRDS(integrated_seurat, paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  } else {
    cat("\nYou missed annotating these clusters:\t", setdiff(list_1, list_2))
    cat("\nThese clusters are not present in data:\t", setdiff(list_2, list_1))
  }
}

#******************************************************************************#
#                         STEP 11: VISUALIZE                                   #
#******************************************************************************#

# Visualize subtypes
for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")){
  split <- NULL
  visualize_UMAP()
}

# Re-annotate full seurat object & visualize
celltypes <- c("Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")
re_annotate(celltypes)
celltype <- NULL
split <- "Condition"
visualize_UMAP()

# Plot dotplots of markers
for (celltype in c("All Markers", "Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")){
  visualize_dotplot()
}

# #****************PLOT STACKED BAR CHARTS OF CELLTYPES PER SAMPLE***************#
# 
# # Create a workbook to store the results
# wb <- openxlsx::createWorkbook()
# 
# for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Endothelial", "Lymphoid")){
#   
#   # Load the integrated seurat object
#   integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
#                                       dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
#   
#   # Remove unwanted cells and samples
#   integrated_seurat <- subset(x = integrated_seurat,
#                               subset = (cell_class %in% c("Mixed", "Unclassified") | (Condition %in% c("Normal"))),
#                               invert = TRUE)
#   
#   # Calculate cells per cluster
#   cells_per_cluster <- integrated_seurat@meta.data %>% 
#     dplyr::group_by(Sample, sub_type) %>%
#     dplyr::count() %>%
#     tidyr::pivot_wider(id_cols = Sample, names_from = sub_type, values_from = n) %>%
#     base::replace(is.na(.), 0)
#   
#   # Calculate percent of cells in each cluster for each sample
#   # Divide by rowSums to adjust for difference in cell number between samples
#   cells_per_cluster_percent <- cells_per_cluster %>%
#     dplyr::mutate(across(.cols = everything())*100/rowSums(across()))
#   
#   # Calculate total cells in each cluster across all samples
#   cells_per_cluster_total <- c(list(Sample = "Total Cells per cluster"), colSums(cells_per_cluster[-1]))
#   
#   # Merge all data
#   cells_per_cluster <- dplyr::bind_rows(cells_per_cluster, data.frame(data = " "), 
#                                         cells_per_cluster_percent, data.frame(data = " "),
#                                         cells_per_cluster_total) %>%
#     dplyr::select(-data)
#   
#   # Plot stacked histogram
#   # Calculate percentages
#   data <- cells_per_cluster_percent %>%
#     tidyr::pivot_longer(cols = !Sample, names_to = "sub_type", values_to = "percent")
#   
#   # Plot stacked bar chart
#   ggplot2::ggplot(data = data, aes(x = Sample, y = percent, fill = sub_type)) +
#     geom_bar(stat = "identity", position = "stack", width = 0.95) +
#     theme_classic() + 
#     scale_fill_manual(values = my_palette) +
#     #geom_text(aes(x=Sample, label=percent), position=position_stack(vjust=0.5), fontface="bold", colour="white", size=6, check_overlap=TRUE) +
#     ggplot2::labs(title = "",
#                   fill = "Subtype",
#                   x = "",
#                   y = "% of Cells") +
#     my_theme
#   
#   # Save the plot
#   ggplot2::ggsave(filename = paste0("Subtype_Percent_", celltype, ".pdf"),
#                   plot = last_plot(),
#                   device = "pdf",
#                   path = seurat_results,
#                   scale = 1,
#                   #width = 4*1+5,
#                   #height = 7,
#                   units = c("in"),
#                   dpi = 600,
#                   limitsize = TRUE,
#                   bg = NULL)
#   
#   # Save data
#   openxlsx::addWorksheet(wb = wb, sheetName = celltype)
#   openxlsx::writeData(wb = wb, sheet = celltype, x = cells_per_cluster)
# }
# openxlsx::saveWorkbook(wb, file = paste0(seurat_results, "Subtype_Percent.xlsx"),
#                        overwrite = TRUE)
# 
# #********************PLOT PIE CHART OF CELLTYPES PER SAMPLE********************#
# 
# celltype <- NULL
# 
# # Load the integrated seurat object
# integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
#                                     dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
# 
# # Remove unwanted cells and samples
# integrated_seurat <- subset(x = integrated_seurat,
#                             subset = (cell_class %in% c("Mixed", "Unclassified") | (Condition %in% c("Normal"))),
#                             invert = TRUE)
# 
# # Create a function to generate pie chart for one sample at a time
# piechart <- function(sample_id){
#   
#   # Calculate % of each cell type for sample being analyzed
#   data <- integrated_seurat@meta.data %>%
#     dplyr::filter(Sample == sample_id) %>%
#     dplyr::count(cell_type) %>%
#     dplyr::mutate(percent = round(100*n/sum(n, na.rm=TRUE), digits = 0), label_percent = paste0(percent,"%")) %>%
#     dplyr::arrange(cell_type)
#   
#   # If you run ggplot without themevoid(), you will see 0 and 100% dont overlap.
#   ggplot(data = data, aes(x = "", y = percent, fill = cell_type)) +
#     geom_bar(stat = "identity", width = 3, color = "white") +
#     coord_polar(theta = "y", start = 0, direction = -1) +
#     geom_text(aes(x = 3.2, label = label_percent), position = position_stack(vjust = 0.5), color = "black", size = 3.5, check_overlap = TRUE) +
#     scale_fill_manual(values = my_palette,
#                       aesthetics = "fill") +
#     ggplot2::labs(title = paste0(sample_id), #, " (", sum(data$n) , ")"),
#                   fill = "Clusters",
#                   x = "",
#                   y = "") +
#     theme_void() +        #remove background, grid, numeric labels
#     my_theme +
#     ggplot2::theme(axis.text.x =  element_blank(),
#                    axis.text.y =  element_blank(),
#                    strip.text.x = element_text(family="sans", face="bold",  colour="black", size=10, hjust = 0.5),
#                    legend.position = "none",
#                    axis.line=element_blank(),
#                    axis.ticks=element_blank())
# }
# 
# ncols <- ceiling(length(levels(as.factor(integrated_seurat@meta.data$Sample)))/2)
# 
# pie_plots <- purrr::map(.x = levels(as.factor(integrated_seurat@meta.data$Sample)), 
#                         .f = piechart) %>%
#   cowplot::plot_grid(plotlist = .,
#                      align = "hv",
#                      axis = "tblr",
#                      ncol = dplyr::if_else(ncols <= 6, ncols, 6),
#                      rel_widths = 1,
#                      rel_heights = 1,
#                      labels = NULL,
#                      label_colour = NULL,
#                      label_x = 0,
#                      label_y = 1,
#                      hjust = -0.5,
#                      vjust = 1.5,
#                      scale = 1,
#                      greedy = TRUE,
#                      byrow = TRUE)
# 
# # Extract the legend from one of the plots
# legend_data <- integrated_seurat@meta.data %>%
#   dplyr::filter(Sample == integrated_seurat@meta.data$Sample[1]) %>%
#   dplyr::count(cell_type) %>%
#   dplyr::mutate(percent = round(100*n/sum(n, na.rm=TRUE), digits = 0), label_percent = paste0(percent,"%")) %>%
#   dplyr::arrange(cell_type)
# 
# legend <- cowplot::get_legend(ggplot(data = legend_data, aes(x = "", y = percent, fill = cell_type)) +
#                                 geom_bar(stat = "identity", width = 3, color = "white") +
#                                 coord_polar(theta = "y", start = 0, direction = -1) +
#                                 scale_fill_manual(values = my_palette,
#                                                   aesthetics = "fill") +
#                                 ggplot2::labs(fill = "") +
#                                 my_theme +
#                                 ggplot2::theme(axis.text.x =  element_blank(),
#                                                axis.text.y =  element_blank(),
#                                                strip.text.x = element_text(family="sans", face="bold",  colour="black", size=10, hjust = 0.5),
#                                                legend.position = "bottom",
#                                                legend.direction = "horizontal",
#                                                axis.line=element_blank(),
#                                                axis.ticks=element_blank()))
# 
# # Add legend to bottom of plots
# pie_plots <- cowplot::plot_grid(pie_plots, legend,
#                                 ncol = 1,
#                                 rel_heights = c(1, .1))
# 
# # Save the plot
# ggplot2::ggsave(filename = "Pie_chart_Percent_Cell_Types.pdf",
#                 plot = pie_plots,
#                 device = "pdf",
#                 path = seurat_results,
#                 scale = 1,
#                 width = 11,       #1.5 inch per pie
#                 height = 9,     #2 inch per row + 1 inch for legend
#                 units = c("in"),
#                 dpi = 300,
#                 limitsize = TRUE,
#                 bg = "white")
# 
# # # Piechart looks good but (i) no legend (ii) unable to place multiple plots
# # # in same page
# # pdf(filename = paste0(results_path, "Percent_Cell_Types.pdf"),
# #      width = 8.5,
# #      height = 11,
# #      bg = "white")
# #
# # pie(x = data$percent,
# #     labels = data$label_percent,
# #     border="white",
# #     col=brewer.pal(12,"Paired"))
# #
# # dev.off()
# 
# ########### Box plot GeoMx
# 
# # Read DEG file for each immune cell type
# epi<- read.xlsx(paste0(seurat_results, "Results_Epithelial_id_Male_vs_Female_DEGs.xlsx"))
# fib <- read.xlsx(paste0(seurat_results, "Results_Fibroblasts_id_Male_vs_Female_DEGs.xlsx"))
# t <- read.xlsx(paste0(seurat_results, "Results_Lymphoid - T_id_Male_vs_Female_DEGs.xlsx"))
# b <- read.xlsx(paste0(seurat_results, "Results_Lymphoid - B_id_Male_vs_Female_DEGs.xlsx"))
# NK <- read.xlsx(paste0(seurat_results, "Results_Lymphoid - NK_id_Male_vs_Female_DEGs.xlsx"))
# mac <- read.xlsx(paste0(seurat_results, "Results_Myeloid - Macrophages, DCs_id_Male_vs_Female_DEGs.xlsx"))
# mdsc <- read.xlsx(paste0(seurat_results, "Results_Myeloid - MDSC_id_Male_vs_Female_DEGs.xlsx"))
# 
# colnames(epi)[1] <- "gene"
# colnames(fib)[1] <- "gene"
# colnames(t)[1] <- "gene"
# colnames(b)[1] <- "gene"
# colnames(NK)[1] <- "gene"
# colnames(mac)[1] <- "gene"
# colnames(mdsc)[1] <- "gene"
# 
# genes <- c("Bcl2l1","Cd27","Ctla4", "Gzmb", "Esr1", "Bad", "Cd163", 
#            "Cd44", "Ifngr1", "Pmel", "Il7r", "Cd4", "Cd14", "Casp3", "Cd74")
# # genes <- c("BCL2L1","CD27","CTLA4", "GZMB", "ESR1", "BAD", "CD163", 
# #            "CD44", "IFNGR1", "PMEL", "IL7R", "CD4", "CD14", "CASP3", "CD74")
# 
# epi <- epi %>% dplyr::filter(gene %in% genes) %>% dplyr::filter(padj < 0.05)
# fib <- fib %>% dplyr::filter(gene %in% genes) %>% dplyr::filter(padj < 0.05)
# t <- t %>% dplyr::filter(gene %in% genes) %>% dplyr::filter(padj < 0.05)
# b <- b %>% dplyr::filter(gene %in% genes) %>% dplyr::filter(padj < 0.05)
# NK <- NK %>% dplyr::filter(gene %in% genes) %>% dplyr::filter(padj < 0.05)
# mac <- mac %>% dplyr::filter(gene %in% genes) %>% dplyr::filter(padj < 0.05)
# mdsc <- mdsc %>% dplyr::filter(gene %in% genes) %>% dplyr::filter(padj < 0.05)
# 
# dim(epi)
# dim(fib)
# dim(t)
# dim(b)
# dim(NK)
# dim(mac)
# dim(mdsc)

#***************(RECOMMENDED): REMOVING BATCH SPECIFIC CLUSTERS**************#

# # Sometimes, you will see sample/batch specific clusters where cells of a 
# # particular sample dominate. Rather than removing all cells from such a cluster
# # we can remove ONLY the cells belonging to dominating sample at all 
# # resolutions and re-run entire pipeline for better clustering.
# 
# # In example below, we can just remove cells from N6-0, N5-1, N5-11 and N4-12
# 
# # Sample  0	            1	            10	          11  	       12
# # N1	    40	          38	          240	          0	            130
# # N2	    20	          222	          230	          4   	        78
# # N3	    632	          350	          2639	        8	            325
# # N4	    13	          34	          872	          0	            "5190"
# # N5	    399	          "28556"       487	          "6241"	      96
# # N6	    "52043"	       562	        1900	        47	          300
# # 
# # N1	    0.536408743	  0.509588306	  3.218452461	  0	            1.743328416
# # N2	    0.109110747	  1.211129296	  1.254773595	  0.021822149	  0.425531915
# # N3	    0.854920528	  0.473452824	  3.569834292	  0.010821779	  0.439634765
# # N4	    0.098776689	  0.258339032 	6.62563635	  0	            39.43469341
# # N5	    0.56627874	  40.52795913	  0.691172296	  8.857507806	  0.136247516
# # N6	    71.15045458	  0.768336865	  2.597580149	  0.06425593	  0.410144234
# # 
# # Total   53147	        29762	        6368	        6300	        6119 
# 
# # Skip round 2 if no bad cells are present
# if (nrow(bad_cells) == 1){
#   break
# } else{
#   # Remove dominant and isolated cells
#   filtered_seurat <- subset(x = integ_data,
#                             Cell %in% bad_cells$Cell, 
#                             invert = TRUE)
#   
#   # Change default assay to RNA and remove other assays.
#   # You need to perform SCTransform() again as removing cells will alter 
#   # Pearson residuals -> UMI corrected value for some genes in some cells that
#   # were earlier 0.
#   DefaultAssay(object = filtered_seurat) <- "RNA"
#   filtered_seurat[["SCT"]] <- NULL   
#   filtered_seurat[["integrated"]] <- NULL
# }