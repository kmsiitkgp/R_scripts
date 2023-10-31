#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

#***************************DEFINE ENSEMBL VARIABLES***************************#

# Choose which database you want to get gene annotations from.
# NOTE: Most RNA Seq data has Ensembl ids as rownames
database <- "Entrez"
database <- "Ensembl"

# Choose if data is from human or mice. We will adjust gene names accordingly.
species <- "Homo sapiens"
species <- "Mus musculus"

#****************************DEFINE DESEQ2 VARIABLES***************************#

# NOTE: Make sure there are no white spaces in Target and Reference columns
# in excel file. R will change any white space (" ") to dot ("."). So, replace 
# white space (" ") with underscore ("_") before importing into R.
# NOTE: Metadata MUST have a column named "Batch" before importing into R.
# Example: If samples were isolated on different days & prepared using different
# kits, "Batch" column must have values like "1_Ribo", "1_Poly-A", "2_Ribo", etc

# Define DESeq2 thresholds
padj.cutoff <- 0.1
lfc.cutoff <- 0     # 0.58 for more strict analysis

# Indicate if you want to plot heatmap, volcano plot, PCA plot
heatmap_plot <- FALSE
volcano_plot <- FALSE
PCA_plot <- FALSE
cor_plot  <- FALSE

# Folder name corresponding to project containing input files
proj <- "TCGA"
proj <- "IMVigor210"
proj <- "Boopati_MB49Sigma_DDR2KD3"
proj <- "Boopati_MB49Y-_DDR2KO"
proj <- "Mukta_NA13_MB49_CDH12OE"
proj <- "Mukta_GSE95097"
proj <- "Hany_Y"
proj <- "Hany_YKO"
proj <- "Hany_KO"
proj <- "Hany_OE"
proj <- "Hany_OE_Tumor"
proj <- "TRAMP_GSE79756"
proj <- "Hany_Y_Tumor"
proj <- "Prince_DepMap"

# Define Target and Reference base don project
if (proj =="Mukta_GSE95097"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("8hr", "4hr"),
                      Reference = c("0hr", "0hr"))
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj == "Boopati_MB49Y-_DDR2KO"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Ddr2KO"),
                      Reference = c("Control"))
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj == "Boopati_MB49Sigma_DDR2KD3"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Ddr2KD"),
                      Reference = c("Control"))
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Mukta_NA13_MB49_CDH12OE"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("OE"),
                      Reference = c("Control"))
  Variable2 <- "Cell_line"
  Variable2_value <- "NA13"  #MB49"
}
if (proj=="Hany_KO"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Kdm5d_KO", "Uty_KO"), 
                      Reference = c("scr_KO", "scr_KO"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_YKO"){
  Variable <- "Condition"
  Comparisons <- list(Target =    c("YKO_centromere"), 
                      Reference = c("scr_KO_centromere"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_Y"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Y_neg"),
                      Reference = c("Y_pos"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_Y_Tumor"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Y_pos_Tumor", "Hany_Tumor"),
                      Reference = c("Y_neg_Tumor", "Sigma_Tumor"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_OE"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Kdm5d_OE", "Uty_OE"),
                      Reference = c("scr_OE", "scr_OE"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_OE_Tumor"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Kdm5d_OE_Tumor", "Uty_OE_Tumor"),
                      Reference = c("scr_OE_Tumor", "scr_OE_Tumor"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="TRAMP_GSE79756"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("PAX8-NFE2L2 fusion"),
                      Reference = c("Empty vector"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="TCGA"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Male"),
                      Reference = c("Female"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Prince_DepMap"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Ypos", "Female"),
                      Reference = c("Yneg", "Yneg"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}

if (proj=="scRNASeq_BBN_Nude"){
  Variable <- "Sex"
  Comparisons <- list(Target = c("Male"),
                      Reference = c("Female"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}

if (proj=="scRNASeq"){
  Variable <- "Ystatus"
  Comparisons <- list(Target = c("Yneg"),
                      Reference = c("Ypos"))    
  Variable2 <- "Condition"
  Variable2_value <- "Tumor"
}

#***************************DEFINE HEATMAP VARIABLES***************************#

# Define if data is log transformed already
already_log <- FALSE

# Define if data is scaled already
already_scaled <- FALSE

# Define if you want to perform unsupervised row and column clustering
# NOTE: If set to FALSE, samples (columns) and genes(rows) will be arranged 
# in alphabetical order (default) in heatmap. If you want to arranged in 
# specific order, define below.
row_clustering <- TRUE    # Usually TRUE
col_clustering <- TRUE    # Usually FALSE

# Define if you want genes or samples to be arranged in alphabetical order
# NOTE: If set to FALSE, write the plot_genes in order you want in heatmap
# NOTE: If row_clustering==TRUE, then row_clustering_alphabetical is irrelevant
row_clustering_alphabetical <- TRUE
col_clustering_alphabetical <- FALSE

# List annotations you want on heatmap
# NOTE: anno_columns MUST match one of the column names in metadata
anno_columns <- c("Condition")

# Define colors for heatmap
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(200)
#my_palette <- viridis(200)

# Dimensions of each box in heatmap
cell_width <- NA
cell_height <- NA

#*************************DEFINE VOLCANO PLOT VARIABLES************************#

# Define the target and reference groups
# NOTE: These values MUST be present in "Condition" column of metadata
Target <- Comparisons$Target
Reference <- Comparisons$Reference

# Define if you want to color by padj values or direction (up vs down)
color_by <- "Significance"
color_by <- "Direction"

# Define cutoffs for padj and log2FC
padj_cutoff <- 0.05
log2_cutoff <- 0.58

#******************************************************************************#
#                            STEP 1: GET ANNOTATIONS                           #
#******************************************************************************#

annotations <- get_annotations(species)

# parent directory : directory where input files, results, etc are stored
parent_path <- paste0("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/", proj, "/")
results_path <- paste0("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/", proj, "/")

#******************************************************************************#
#                           STEP 2: IMPORT META DATA                           #
#******************************************************************************#

# NOTE: You have to prepare a xlsx file where 
# column 1 has filenames of HTseq output files, 
# column 2 has Condition and rest of columns with other information.

# Import meta data and add sample names as rownames to meta_data
meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Metadata.xlsx"))

#******************************************************************************#
#                  STEP 3A: IMPORT READ DATA FROM SINGLE FILE                  #                                           
#******************************************************************************#

if (proj %in% c("TCGA", "Prince_DepMap")){
  # Import counts, remove gene names with NA values, add gene names as rownames
  read_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Readdata.xlsx"))
}

#******************************************************************************#
#                 STEP 3B: IMPORT READ DATA FROM MULTIPLE FILES                #                                           
#******************************************************************************#

if (proj != "TCGA"){
  # Create a list of all txt files within folder that will be analyzed
  count_folder <- paste0(parent_path, "HTSeq_count_results/")
  files <- list.files(path=count_folder)
  #files <- intersect(files, rownames(meta_data %>% dplyr::filter(id != "ignore")))
  
  # Create an empty dataframe with 0's
  read_data <- data.frame(0)
  
  # Create the reads table 
  for (i in 1:length(files)){
    
    # Read the txt file
    temp_file <- read.table(file = paste0(count_folder, files[i]), header = FALSE, sep = "\t")     
    
    # The 1st Column will have Ensembl ids. 
    # Gene counts may be in 2nd or 3rd column. Append appropriate column.
    if (sum(temp_file[2], na.rm = TRUE) == 0 & sum(temp_file[3], na.rm = TRUE) > 0){
      read_data <- bind_cols(read_data, temp_file[,3])
    } else if (sum(temp_file[2], na.rm = TRUE) > 0 & sum(temp_file[3], na.rm = TRUE) ==0){
      read_data <- bind_cols(read_data, temp_file[,2])                  
    } else{
      print("Error: Gene counts NOT PRESENT in either column 2 or 3 of count file")
    }
    
    # Rename the column names to sample names
    colnames(read_data)[i+1] <- gsub(pattern="\\..*$", replacement="", x=files[i])
  }
  
  # Check if all count files have same order of genes in the rows so that files can be merged together
  gene_list <- temp_file[, 1]
  for (i in 1:length(files)){
    
    temp_file <- read.table(file = paste0(count_folder, files[i]), header = FALSE, sep = "\t")
    genes <- temp_file[, 1]
    
    if (identical(gene_list, genes)){ 
      print("Gene order is same between count files")
    } else {
      print("Gene order is different between the count files")
    }
  }
  
  # Add gene names to 1st column
  read_data[, 1] <- gene_list
  colnames(read_data)[1] <- "SYMBOL"
  
  # Last 5 rows in HTSeq count output are __no_feature,  __ambiguous, 
  # __too_low_aQual, __not_aligned, __alignment_not_unique. So, we remove them.
  read_data <- read_data[1:(nrow(read_data)-5),]
}

analyze_DESeq2()

# #*******************************DIAGNOSTIC TESTS*******************************#
# 
# celltype <- NULL
# 
# # (i) To view counts of specific gene across samples
# plotCounts(dds, gene=which.min(res$padj), intgroup=Variable)           # gene with lowest padj
# plotCounts(dds, gene=which.min(res$log2FoldChange), intgroup=Variable) # gene with lowest log2FC
# plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup=Variable) # gene with highest log2FC
# plotCounts(dds, gene="ENSMUSG00000030598", intgroup=Variable)          # a specific gene
# 
# # (ii) Plot MA
# # The output of plotMA is not a ggplot2 object. So, we cant use ggsave()
# pdf(filename = paste0(diagnostics_path, "Diagnostic_MA_Plot_", celltype, ".pdf"),
#     width = 8.5,
#     height = 11,
#     units = "in",
#     quality = 75,
#     bg = "white",
#     res = 600)
# 
# # Blue points indicate genes with adj p value <0.1
# plotMA(object = res,
#        alpha = 0.1,
#        main = "",
#        xlab = "Mean of Normalized Counts",
#        #ylim = c(-5,5), 
#        MLE = FALSE) 
# 
# dev.off()
# # To identify the genes interactively, run the 2 lines below. 
# # Then click on multiple dots and click Finish. A list of genes will be displayed 
# # idx <- identify(res$baseMean, res$log2FoldChange)
# # rownames(res)[idx]
# 
# # (iii) Plot PCA using rld or vst
# # Variable is defined at beginning of DESeq2. So, we set intgroup=Variable.
# rld <- rlog(object = dds,
#             blind = TRUE,
#             # intercept = ,
#             # betaPriorVar = ,
#             fitType = "parametric")
# 
# vst <- varianceStabilizingTransformation(object = dds,
#                                          blind = TRUE,
#                                          fitType = "parametric")
# 
# for (process in c("rld", "vst")){
#   
#   plotPCA(object = get(process),
#           intgroup = Variable,
#           ntop = 500,
#           returnData = FALSE) +
#     geom_text(aes(label = name), nudge_x = 2, nudge_y = 2)
#   
#   ggplot2::ggsave(filename = paste0("Diagnostic_PCA_plot_using_", process, "_", celltype, ".pdf"),
#                   plot = last_plot(),
#                   device = "pdf",
#                   path = diagnostics_path,
#                   scale = 1,
#                   #width = 8.5,
#                   #height = 11,
#                   units = c("in"),
#                   dpi = 300,
#                   limitsize = TRUE,
#                   bg = "white")
# }
# 
# # (iv) Hierarchical clustering of samples using rld or vst
# rld_mat <- DESeq2::assay(rld)     # extract the matrix from a DESeq2 object
# rld_cor <- cor(x = rld_mat,       # compute pairwise correlation values
#                y = NULL,
#                use = "everything",
#                method = "pearson") 
# 
# # Check the output of cor(), make note of the rownames and colnames
# head(rld_cor)
# 
# vst_mat <- assay(vst)  
# vst_cor <- cor(vst_mat) 
# 
# # Check the output of cor(), make note of the rownames and colnames
# head(vst_cor)    
# 
# for (process in c("rld", "vst")){
#   
#   pheatmap::pheatmap(mat = get(paste0(process, "_cor")),
#                      color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100),
#                      breaks = NA, 
#                      border_color = "white", #"grey60",
#                      cellwidth = NA, 
#                      cellheight = NA, 
#                      scale = "none",   
#                      cluster_rows = TRUE,   #cluster the rows
#                      cluster_cols = TRUE,   #cluster the columns
#                      clustering_distance_rows = "euclidean",
#                      clustering_distance_cols = "euclidean",
#                      clustering_method = "complete",
#                      legend = TRUE, 
#                      legend_breaks = NA,
#                      legend_labels = NA, 
#                      #annotation_row = ,  
#                      #annotation_col = , 
#                      annotation_colors = dplyr::if_else(nrow(col_annotation)+nrow(row_annotation) > 0, ann_colors, NA),
#                      annotation_legend = TRUE,
#                      annotation_names_row = TRUE,
#                      annotation_names_col = TRUE,
#                      show_rownames = dplyr::if_else(nrow(mat)<80, TRUE, FALSE, missing = NULL), 
#                      show_colnames = dplyr::if_else(ncol(mat)<30, TRUE, FALSE, missing = NULL),
#                      fontsize = 8, 
#                      fontsize_row = 8, 
#                      fontsize_col = 8,
#                      angle_col = c("270", "0", "45", "90", "315"),
#                      fontsize_number = 0.8*fontsize, 
#                      #labels_row = display_row,
#                      #labels_col = display_col,
#                      filename = paste0(diagnostics_path, "Diagnostic_Correlation_Heatmap_using_", process, "_", celltype, ".pdf"))
# }
# 
# # (v) Checking if mean < variance (for NB model) or Mean = Variance (for Poisson model). 
# # Each point is a gene denoted by (x,y) where 
# # x = mean_count of gene across all samples & y = variance_count of gene across all samples
# 
# mean_counts <- apply(read_data[,], 1, mean)
# variance_counts <- apply(read_data[,], 1, var)
# df <- data.frame(mean_counts, variance_counts)
# ggplot(df) +
#   geom_point(aes(x=mean_counts, y=variance_counts)) +
#   geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
#   scale_y_log10() +
#   scale_x_log10()
# 
# ggplot2::ggsave(filename = paste0("Diagnostic_Scatter_plot_", celltype, ".pdf"),
#                 plot = last_plot(),
#                 device = "pdf",
#                 path = diagnostics_path,
#                 scale = 1,
#                 #width = 8.5,
#                 #height = 11,
#                 units = c("in"),
#                 dpi = 300,
#                 limitsize = TRUE,
#                 bg = "white")
# 
# # (vi) Plot dispersion estimates. Higher the mean, lower the dispersion
# # The output of plotDispEsts() is not a ggplot2 object. So, we cant use ggsave()
# 
# pdf(filename = paste0(diagnostics_path, "Diagnostic_Dispersion_Estimates_", celltype, ".pdf"),
#     width = 8.5,
#     height = 11,
#     units = "in",
#     quality = 75,
#     bg = "white",
#     res = 600)
# 
# plotDispEsts(object = dds,
#              #ymin,
#              genecol = "black",
#              fitcol = "red",
#              finalcol = "dodgerblue",
#              legend = TRUE,
#              xlab = "Mean of Normalized counts" ,
#              ylab = "Dispersion",
#              log = "xy",
#              cex = 0.45)
# 
# dev.off()
# 
# # # (vii) Plot p-value histogram
# # hist(res$pvalue, col="lightblue", breaks=20)
# # 

# # 
# # # (ix) Plot a histogram for one of the samples to see how the counts are distributed. Adjust xlim if needed
# # ggplot(read_data, aes(x = results.S10_R1_trimmed.fastq.gz.csv)) +
# #   geom_histogram(stat = "bin", bins = 200) +
# #   xlim(-5,1000) +
# #   xlab("Raw expression counts") +
# #   ylab("Number of genes") +
# #   scale_color_brewer(palette="Dark2")
# # 
# # # (x) Extract counts, size factors, etc.. from dds
# # dds <- estimateSizeFactors(dds)         #redundant if you already ran dds <- DESeq(dds)
# # counts <- counts(dds)                   #counts[i,j] = raw_count of gene i in sample j
# # sizefactors <- sizeFactors(dds)         #sizefactors[j] = median (counts[,j]/geometric_mean(counts[i,]))
# # colSums(counts(dds))                    #Total number of raw counts per sample
# # colSums(counts(dds, normalized=TRUE))   #Total number of normalized counts per sample