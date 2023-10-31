#!/usr/bin/env Rscript

# NOTE: Some packages like Seurat-disk, ScopeLoomR need hdf5r package which in 
# turn needs HDF5 library files to install as well as even work. Similarly, 
# metap needs fftw3. So, 
# (i) qrsh -> conda activate R -> module load hdf5/1.8.18 -> R -> install... 
# (ii) qrsh -> conda activate R -> conda install -c conda-forge fftw -> R -> install...

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
library("TCGAbiolinks")         # Needed for TCGA data analysis
library("ensembldb")            # Needed for annotating genes
library("AnnotationHub")        # Needed for annotating genes
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("affy")                 # Needed for micro-array analysis
library("lumi")                 # Needed for micro-array analysis
#library("illuminaHumanv4.db")  # Needed for micro-array analysis
library("limma")                # Needed for micro-array analysis
library("ChIPQC")               # Needed for ChIP analysis
library("fgsea")                # Needed for GSEA analysis
library("enrichplot")           # Needed for GSEA analysis
library("clusterProfiler")      # Needed for GSEA analysis
#library("pathview")            # Needed for GSEA analysis
library("msigdbr")              # Needed for GSEA analysis
library("DESeq2")               # Needed for Differential Expression analysis
library("progeny")
library("dorothea")
library("viper")
library("infercnv")

# Data wrangling packages
library("openxlsx")             # Needed for reading, writing xlsx files
library("dplyr")                # Needed for data wrangling
library("tibble")               # Needed for data wrangling
library("stringr")              # Needed for data wrangling
library("purrr")                # Needed for data wrangling

# Graph plotting packages
library("ggplot2")              # Needed for making graphs
library("cowplot")              # Needed for merging multiple graphs
library("viridis")              # Needed for nice graph coloring
library("RColorBrewer")         # Needed for nice graph coloring
library("ggrepel")              # Needed for making graphs prettier
library("ggpubr")               # Needed for adding p values to graphs
library("ggbeeswarm")           # Needed for proper positioning of labels in scatter plots
library("colorspace")

# Specialized Graph plotting packages
library("pheatmap")             # Needed for making heatmap
library("ggridges")             # Needed for making ridgeplots
#library("hrbrthemes")           # Needed for modern look of plots
library("VennDiagram")          # Needed for making Venn diagram
library("survival")             # Needed for making survival curves
library("survminer")            # Needed for making survival curves

# Single cell analysis packages
library("Seurat")               # Needed for single cell analysis
library("SeuratDisk")           # Needed for reading h5ad files
library("SCopeLoomR")           # Needed for reading loom files
library("SCENIC")               # Needed for SCENIC analysis
library("DoubletFinder")        # Needed for identifying doublets
library("Augur")
library("ktplots")              # Needed for plotting cellphonedb results
library("CellChat")
library("patchwork")
library("wordcloud")

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

# Choose project specific folder where data will be saved
proj <- "scRNASeq_Chen" 
# "scRNASeq_BBN_C57B6"      
# "scRNASeq_BBN_Rag"
# "scRNASeq_NA13_CDH12_C57B6"
# "scRNASeq_NA13_CDH12_Nude"
# "scRNASeq_GSE164557"
# "scRNASeq_GSE132042"
# "scRNASeq_Chen"
# "scRNASeq_Simon"

# Choose if data is from human or mice. We will adjust gene names accordingly.
species <- "Homo sapiens"
species <- "Mus musculus"

# Choose if current data is h5ad file or output of cellranger count
h5ad_file <- "yes"
h5ad_file <- "no"

# Choose if current data was demultiplexed.
demultiplexed <- "yes"
demultiplexed <- "no"

# Choose if you want to generate whitelist of barcodes for CITE-Seq-count
whitelist <- "yes"
whitelist <- "no"

# Choose which database you want to get gene annotations from
database <- "Entrez"
database <- "Ensembl"

# Choose xlsx file with metadata for the project. 
# NOTE: This must be an xlsx file with one column as "Sample" which matches with
# Seurat metadata column "Sample" so that we can merge additional info from 
# this file with Seurat object. It MUST also be named using same value you
# set to "proj" variable above and a suffix of "_Metadata". 
metafile <- paste0(proj, "_Metadata.xlsx")

# Choose reference samples to be used during integration. 
# NOTE: Use batch numbers as references, not individual patients.
# Set the column in meta.data and its value to be used as reference
# For example, to use all samples having value as "Female" under column "Sex" &
# "BBN" under column "Treatment" as reference, set 
# ref1 <- "Sex"; ref1_value <- "Female"; ref2 <- "Treatment"; ref2_value <- "BBN"
ref1 <- NA            
ref1_value <- NA         
ref2 <-  NA            
ref2_value <- NA           

# Choose the variables to use in DESeq2
padj.cutoff <- 0.1
lfc.cutoff <- 0       # set to 0.58 for more strict analysis
Variable <- "Sex"
Target <- "Male"
Reference <- "Female"

# Choose the variable to subset. 
# For example, if the seurat object has "BBN" treated samples and "control" 
# samples but you want to analyze ONLY BBN samples, set those variables here.
Variable2 <- "Condition"
Variable2_value <- "Tumor"

# Change default limit for allowable object sizes within R 
options(future.globals.maxSize = 30000 * 1024^2)

# scripts directory        : contains scripts, xlsx files, etc..
# parent directory         : contains all sub-directories
# feature matrix directory : contains raw feature matrices of cellranger count
# hto matrix directory     : contains hto umicount output of CITE-Seq-count
# seurat directory         : contains singlet seurat objects after hashtag demultiplexing
# diagnostics directory    : contains QC and diagnostic plots
# results directory        : contains intermediate rds files, final figures

# NOTE: The matrix directory should have the raw_feature_bc_matrix folders for
# each sample. Each of these raw_feaure_bc_matrix folders SHOULD have 
# matrix.mtx, genes.tsv (or features.tsv) & barcodes.tsv files.
# You need to manually add sample name as prefix to each of the
# raw_feature_bc_matrix folders before starting analysis.
# Eg: FB1_raw_feature_bc_matrix, FB2_raw_feature_bc_matrix,.....

scripts_path <-         "/hpc/home/kailasamms/projects/scRNASeq/"
parent_path <-          paste0("/hpc/home/kailasamms/scratch/", proj, "/")
feature_matrix_path <-  paste0(parent_path, "raw_feature_bc_matrix/")
hto_matrix_path <-      paste0(parent_path, "raw_hto_bc_matrix/")
demux_results <-        paste0(parent_path, "results_demux/")
seurat_results <-       paste0(parent_path, "results_seurat/")
pyscenic_results <-     paste0(parent_path, "results_pyscenic/")
scvelo_results <-       paste0(parent_path, "results_scvelo/")
velocyto_results <-     paste0(parent_path, "results_velocyto/")
cellphonedb_results <-  paste0(parent_path, "results_cellphonedb/")
cellchat_results <-     paste0(parent_path, "results_cellchat/")
dorothea_results <-     paste0(parent_path, "results_dorothea/")
diagnostics_path <-     paste0(parent_path, "diagnostics/")

# Load the cell cycle markers
cell_cycle_markers <- openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, "Cell_Cycle_Markers.xlsx"),
                                          colNames = TRUE)
# Create a list of S and G2M markers
if (species == "Homo sapiens"){
  cell_cycle_genes <- cell_cycle_markers$Human_Gene
  s_genes <- cell_cycle_markers$Human_Gene[which(cell_cycle_markers$Phase=="S")]
  g2m_genes <- cell_cycle_markers$Human_Gene[which(cell_cycle_markers$Phase=="G2/M")]
} else {
  cell_cycle_genes <- cell_cycle_markers$Mouse_Gene
  s_genes <- cell_cycle_markers$Mouse_Gene[which(cell_cycle_markers$Phase=="S")]
  g2m_genes <- cell_cycle_markers$Mouse_Gene[which(cell_cycle_markers$Phase=="G2/M")]
}

# Define aspect ratio
aspect_ratio <- 1

# Define axis font etc to use in all plots
my_theme <- ggplot2::theme(aspect.ratio = aspect_ratio,
                           plot.title =   element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
                           axis.title.x = element_text(family="sans", face="bold",  colour="black", size=12, hjust = 0.5, vjust = 0,   angle = 0),
                           axis.title.y = element_text(family="sans", face="bold",  colour="black", size=12, hjust = 0.5, vjust = 1,   angle = 90),
                           legend.title = element_text(family="sans", face="bold",  colour="black", size=12, hjust = 0,   vjust = 1,   angle = 0),
                           axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5, vjust = 0.5, angle = 45),
                           axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5, vjust = 0.5, angle = 0),
                           legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                           strip.text.x = element_text(family="sans", face="bold",  colour="black", size=10, hjust = 0.5),
                           #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                           legend.position = "right",
                           legend.justification = "left",
                           legend.direction = "vertical",
                           legend.key.height = unit(0.5, 'cm'),
                           legend.key.width  = unit(0.5, 'cm'), 
                           legend.text.align = 0)

# Assign colors for UMAP. The current palette supports up to 33 cell types
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#8C564B", 
                "#E377C2", "#BCBD22", "#17BECF", "#FFC61E", "#762A83",
                "#333333", "#FF1F5B", "#B8E80C")

my_palette <- c(my_palette, 
                colorspace::adjust_transparency(col = my_palette, alpha = 0.3), 
                colorspace::adjust_transparency(col = my_palette, alpha = 0.6))

#******************************************************************************#
#                               GENE ANNOTATIONS                               #
#******************************************************************************#

# This function returns a dataframe with following columns:
# "ENSEMBL_ID", "ENTREZ_ID", "SYMBOL", "SYMBOL_ENTREZ", "BIOTYPE", 
# "BIOTYPE_ENTREZ", "START", "END", "CHR", "STRAND", "DESCRIPTION"

get_annotations <- function(species){
  
  # AnnotationHub has SYMBOL-ENSEMBL_ID info ONLY.
  # AnnotationDbi has SYMBOL-ENSEMBL_ID as well as SYMBOL-ENTREZ_ID info.
  
  #**************************GET ENSEMBL ANNOTATIONS***************************#
  # Connect to AnnotationHub
  ah <- AnnotationHub::AnnotationHub()
  
  # Access the Ensembl database for organism
  ahDb <- AnnotationHub::query(x = ah,
                               pattern = c(species, "EnsDb"),
                               ignore.case = TRUE)
  
  # Acquire the latest annotation files
  id <- ahDb %>%
    mcols() %>%
    rownames() %>%
    tail(n = 1)
  
  # Download the appropriate Ensembldb database
  edb <- ah[[id]]
  
  # Extract gene-level information from database
  ensembl <- ensembldb::genes(x = edb,
                              return.type = "data.frame")
  
  
  # Select annotations of interest
  ensembl <- ensembl %>%
    dplyr::rename(ENSEMBL_ID=gene_id, SYMBOL=gene_name, 
                  BIOTYPE=gene_biotype, START=gene_seq_start, END=gene_seq_end, 
                  CHR=seq_name, STRAND=seq_strand,  DESCRIPTION=description) %>%
    dplyr::mutate(SYMBOL = dplyr::case_when(nchar(SYMBOL) == 0 ~ NA,
                                            TRUE ~ SYMBOL)) %>%
    dplyr::select(ENSEMBL_ID, SYMBOL, BIOTYPE, START, END, CHR, STRAND, DESCRIPTION)
  
  #***************************GET ENTREZ ANNOTATIONS***************************# 
  
  # NOTE: mapIds can ONLY retrieve one of "EMSEMBL/SYMBOL/GENETYPE" at a time
  # mapping <- AnnotationDbi::mapIds(x = org.Hs.eg.db, 
  #                                  keys = keys(org.Hs.eg.db),
  #                                  keytype = "ENTREZID", 
  #                                  column = "SYMBOL") %>%
  #   as.data.frame(do.call(cbind, list(.))) %>%
  #   tibble::rownames_to_column("ENTREZID") %>%
  #   dplyr::rename(ENTREZID = identity(1), SYMBOL = identity(2))
  
  if (species == "Homo sapiens"){
    entrez <- AnnotationDbi::select(x = org.Hs.eg.db, keys = keys(org.Hs.eg.db),
                                    columns = c("ENSEMBL", "SYMBOL","GENETYPE"))
  } else if (species == "Mus musculus"){
    entrez <- AnnotationDbi::select(x = org.Mm.eg.db, keys = keys(org.Mm.eg.db),
                                    columns = c("ENSEMBL", "SYMBOL","GENETYPE"))
  }
  
  colnames(entrez) <- c("ENTREZ_ID", "ENSEMBL_ID", "SYMBOL_ENTREZ", "BIOTYPE_ENTREZ")
  
  
  # Merge ensembl and entrez dataframes
  annotations <- dplyr::full_join(ensembl, entrez, by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
    dplyr::select(ENSEMBL_ID, ENTREZ_ID, SYMBOL, SYMBOL_ENTREZ, BIOTYPE, 
                  BIOTYPE_ENTREZ, START, END, CHR, STRAND, DESCRIPTION)
  
  #DBI::dbDisconnect(conn = ah)
  return(annotations)
}

#******************************************************************************#
#        NORMALIZE DATA, IDENTIFY HIGHLY VARIABLE FEATURES, SCALE DATA,        # 
#                PERFORM DIMENSIONAL REDUCTION USING PCA & UMAP                #
#******************************************************************************#

# Use the sctransform method as a more accurate method of normalizing, 
# estimating the variance of the filtered data, and identifying the most 
# variable genes. By default, sctransform accounts for cellular sequencing 
# depth (i.e. nUMIs). Also, we can regress out variation from cell cycle genes
# and mitochondrial genes if needed. 
# NOTE: We MUST perform SCTransform on each sample INDIVIDUALLY, not on merged
# seurat object.

# Refer https://satijalab.org/seurat/articles/sctransform_vignette.html
# The residuals (normalized values) are stored in pbmc[["SCT"]]@scale.data and 
# used directly as input to PCA. Please note that this matrix is non-sparse, and
# can therefore take up a lot of memory if stored for all genes. To save memory,
# we store these values only for variable genes, by setting the 
# return.only.var.genes = TRUE by default in the SCTransform().

# To assist with visualization and interpretation, we also convert Pearson 
# residuals back to ‘corrected’ UMI counts. You can interpret these as the UMI 
# counts we would expect to observe if all cells were sequenced to the same depth.
# The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@counts. 

# The log-normalized versions of these corrected counts are stored in 
# pbmc[["SCT"]]@data, which are very helpful for visualization.

# You can use the corrected log-normalized counts for differential expression
# and integration. However, in principle, it would be most optimal to perform
# these calculations directly on the residuals (stored in the scale.data slot) 
# themselves. This is not currently supported in Seurat v3, but will be soon.

# NOTE: If you use lapply, even if you return(split_seurat), it wont work 
# because lapply saves S4 objects in base namespace, not within user defined 
# function namespace. So, split_seurat can be accessed as a global variable.

sctransform_data <- function(filtered_seurat){
  
  # Split the merged object into individual objects based on samples
  split_seurat <- Seurat::SplitObject(object = filtered_seurat,
                                      split.by = "Sample")
  
  # Identify samples with less than 50 cells
  samples_with_few_cells <- filt_data@meta.data %>% 
    dplyr::count(Sample) %>% dplyr::filter(n<50)
  
  # Remove samples with less than 50 cells so that RunPCA() doesnt
  split_seurat <- split_seurat[!(names(split_seurat) %in% samples_with_few_cells$Sample)]
  
  # Normalize, identify variable features, SCTransform each dataset independently
  for (i in 1:length(split_seurat)){
    
    # Normalize the data
    split_seurat[[i]] <- Seurat::NormalizeData(object = split_seurat[[i]],
                                               assay = "RNA",
                                               normalization.method = "LogNormalize",
                                               scale.factor = 10000,
                                               margin = 1,
                                               verbose = TRUE)
    
    # SCTransform() is better than FindVariableFeatures() & ScaleData()
    # split_seurat[[i]] <- Seurat::FindVariableFeatures(object = split_seurat[[i]],
    #                                                   assay = "RNA",
    #                                                   selection.method = "vst",
    #                                                   nfeatures = 2000)
    # split_seurat[[i]] <- Seurat::ScaleData(object = split_seurat[[i]],
    #                                        features = NULL,
    #                                        assay = "RNA",
    #                                        vars.to.regress = NULL)
    
    # Perform cell cycle scoring
    split_seurat[[i]] <- Seurat::CellCycleScoring(object = split_seurat[[i]],
                                                  s.features = s_genes,
                                                  g2m.features = g2m_genes,
                                                  ctrl = NULL)
    
    split_seurat[[i]]$CC.Score <- split_seurat[[i]]$G2M.Score-split_seurat[[i]]$S.Score
    
    # Perform scaling & variable feature identification usign SCTransform()
    split_seurat[[i]] <- Seurat::SCTransform(object = split_seurat[[i]],
                                             assay = "RNA",
                                             new.assay.name = "SCT",
                                             reference.SCT.model = NULL,
                                             do.correct.umi = TRUE,
                                             ncells = 5000,
                                             residual.features = NULL,
                                             variable.features.n = 3000,
                                             vars.to.regress = c("CC.Score","MitoRatio"),
                                             do.scale = FALSE,
                                             do.center = TRUE,
                                             return.only.var.genes = TRUE)
    
    # Remove ribosomal, Riken, predicted and mitochondrial genes from
    # VariableFeatures so that PCA, UMAP and hence clustering are not affected
    var_f <- split_seurat[[i]]@assays$SCT@var.features
    var_f <- var_f[!str_detect(string = var_f,
                               pattern = dplyr::if_else(species=="Homo sapiens",
                                                        "^RP[SL]|RIK$|^MT-|^GM[0-9.]+$",
                                                        "^Rp[sl]|Rik$|^mt-|^Gm[0-9.]+$"))]
    split_seurat[[i]]@assays$SCT@var.features <- var_f
    cat("\nFinal number of variable features:", length(var_f), "\n")
    
    # Perform dimensional reduction using PCA on SCT assay variable features
    split_seurat[[i]] <- Seurat::RunPCA(object = split_seurat[[i]],
                                        assay = "SCT",
                                        features = NULL,
                                        ndims.print = 1,
                                        nfeatures.print = 1,
                                        reduction.name = "pca",
                                        reduction.key = "PC_")
    
    # Perform dimensional reduction using UMAP on PCA dimensions
    split_seurat[[i]] <- Seurat::RunUMAP(object = split_seurat[[i]],
                                         dims = 1:40,
                                         reduction = "pca",
                                         reduction.name = "umap",
                                         reduction.key = "UMAP_")
  }
  return(split_seurat)
}

doublet_removal <- function(){

  # Perform doublet removal if needed

  # # Adjust doublet rate based on cell counts in each batch. Higher the doublet
  # # rate more cells will be marked as doublets. B8A and B8B have 30,000+ nuclei
  # # while only 20,000 nuclei were only expected to be recovered. So, using
  # # doublet rate of 0.4, we expect 12000 cells will be marked as doublets.
  # if (i == 8 | i == 9){
  # d_rate = 0.50
  # } else if (i == 1){
  # d_rate = 0.40
  # } else if (i == 4 | i == 5 | i == 6){
  # d_rate = 0.20
  # } else{
  # d_rate = 0.05
  # }

  # cat("\ndoublet rate:", d_rate, "\n")
  # nExp <- round(d_rate*nrow(split_seurat[[i]]@meta.data))

  # ## pK Identification (varies based on data set)
  # sweep.res <- paramSweep_v3(split_seurat[[i]], PCs = 1:40, sct = TRUE)
  # sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  # bcmvn <- find.pK(sweep.stats)
  # print(bcmvn)
  # # bcmvn has factors. So, use as.numeric(as.character())
  # pK =  as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric==max(bcmvn$BCmetric)]))
  # cat("\npK:", pK, "\n")

  # split_seurat[[i]] <- doubletFinder_v3(split_seurat[[i]],
  # PCs = 1:40,
  # pN = 0.25,
  # pK = pK,
  # nExp = nExp,
  # reuse.pANN = FALSE,
  # sct = TRUE)

  # DF.name = colnames(split_seurat[[i]]@meta.data)[grepl("DF.classification", colnames(split_seurat[[i]]@meta.data))]
  # split_seurat[[i]]@meta.data <- split_seurat[[i]]@meta.data %>% rename(DF.classification = DF.name)

  # print(head(split_seurat[[i]]@meta.data %>% count(DF.classification)))

  # DimPlot(split_seurat[[i]],
  # reduction = "umap",
  # group.by = "DF.classification")

  # ggsave(filename = paste0(split_seurat[[i]]@meta.data$Sample[1], "_UMAP_doublet_finder.jpg"),
  # path = seurat_results,
  # bg = "white")

  # split_seurat[[i]] <- subset(x = split_seurat[[i]],
  # DF.classification == "Singlet")

  # print(head(split_seurat[[i]]@meta.data %>% count(DF.classification)))
}

#******************************************************************************#
#                       VISUALIZE DATA BEFORE INTEGRATION                      #
#******************************************************************************#

plot_pre_integration <- function(split_seurat){
  
  # Create a function to plot PCA
  PCA_plot <- function(i){
    Seurat::DimPlot(object = split_seurat[[i]],
                    reduction = "pca",
                    split.by = NULL,
                    label = FALSE,
                    raster = FALSE,
                    ncol = NULL,
                    combine = TRUE) +
      NoLegend() +
      ggplot2::labs(title = names(split_seurat)[i])
  }
  
  # Create a function to plot UMAP
  UMAP_plot <- function(i){
    Seurat::DimPlot(object = split_seurat[[i]],
                    reduction = "umap",
                    split.by = NULL,
                    label = FALSE,
                    raster = FALSE,
                    ncol = NULL,
                    combine = TRUE) +
      NoLegend() +
      ggplot2::labs(title = names(split_seurat)[i])
  }
  
  # Plot PCA & UMAP before integration
  funcs <- c("PCA_plot", "UMAP_plot")
  for (j in 1:length(funcs)){
    
    purrr::map(.x = 1:length(split_seurat), 
               .f = get(funcs[j])) %>%
      cowplot::plot_grid(plotlist = .,
                         align = "hv",
                         axis = "tblr",
                         nrow = NULL,
                         ncol = dplyr::if_else(length(split_seurat) < 3, 2, 3),
                         rel_widths = 1,
                         rel_heights = 1,
                         greedy = TRUE,
                         byrow = TRUE)
    
    # Save the plot
    ggplot2::ggsave(filename = paste0(funcs[j], "_Pre_Integration_", celltype, ".jpg"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = diagnostics_path,
                    scale = 1,
                    width = 8.5,
                    height = 11,
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = "white")
  }
} 

#******************************************************************************#
#               PREPARE THE DATA FOR INTEGRATION & INTEGRATE DATA              #
#******************************************************************************#

# As you see from the UMAP, the cells cluster differently in each sample. 
# To find the same cell population (say macrophages) between 2 samples,
# it is necessary for both samples to have similar clustering pattern in UMAP.
# So, we have to integrate the samples. 

# The goal of integration is to ensure that cell types of one 
# condition/dataset align with the same cell types of the other 
# conditions/datasets (e.g. macrophages in control samples align with 
# macrophages in stimulated condition).

# To integrate, we will use the shared highly variable genes from each 
# condition identified using SCTransform, then, we will "integrate" or 
# "harmonize" the conditions to overlay cells that are similar or have a 
# "common set of biological features" between groups. 

#*********STEP 7A: DECLARE REFERENCE SAMPLES FOR INTEGRATING THE DATA**********#

#***STEP 7B: SELECT 3000 MOST VARIABLE GENES TO USE FOR INTEGRATING THE DATA***#

#******************STEP 7C: FIND RESIDUALS FOR MISSING GENES*******************#
# Each sample has different 3000 most variable genes. Gene X which is most 
# variable among cells of sample A may not be one of the top 3000 most variable 
# genes in sample B. PrepSCTIntegration() will calculate Pearson residuals for
# missing genes so that all samples have the same 3000 genes

#******STEP 7D: FIND COMMON ANCHORS BETWEEN SAMPLES TO INTEGRATE THE DATA******#
# NOTE: Data must be scaled & PCA must have been run before doing cca or rpca
# in this step.

# This step is computationally intensive and involves the following steps:

# (i) Perform canonical correlation analysis (CCA):
# CCA identifies shared sources of variation between the conditions/groups. It
# is a form of PCA, in that it identifies the greatest sources of variation in
# the data, but only if it is shared or conserved across the conditions/groups
# (using the 3000 most variant genes from each sample). This step roughly aligns
# the cells using the greatest shared sources of variation.

# NOTE: The shared highly variable genes are used because they are the most 
# likely to represent those genes distinguishing the different cell types 
# present.

# (ii) Identify anchors or mutual nearest neighbors (MNNs) across datasets 
# (sometimes incorrect anchors are identified): MNNs can be thought of as 
# 'best buddies'. For each cell in one condition:   
# (a) The cell's closest neighbor in the other condition is identified based on
# gene expression values - it's 'best buddy'.
# (b) The reciprocal analysis is performed, and if the two cells are 'best 
# buddies' in both directions, then those cells will be marked as anchors to 
# 'anchor' the two datasets together.

# NOTE: The difference in expression values between cells in an MNN pair 
# provides an estimate of the batch effect, which is made more precise by 
# averaging across many such pairs. A correction vector is obtained and applied
# to the expression values to perform batch correction."
# 
# (iii) Filter anchors to remove incorrect anchors:
# Assess the similarity between anchor pairs by the overlap in their local 
# neighborhoods (incorrect anchors will have low scores)

# NOTE: cca is computationally intensive if more than 2 samples are integrated.
# In such cases use "rpca". Also, using reference based integration is faster.

#******STEP 7E: FIND OPTIMUM k.weight FOR USE IN Seurat::IntegrateData()*****#
# k.weight MUST be less than number of anchors. Else, error will be thrown.

#************************STEP 8A: INTEGRATE THE DATA*************************#
# Use anchors and corresponding scores to transform the cell expression values,
# allowing for the integration of the conditions/datasets (different samples, 
# conditions, datasets, modalities)

# NOTE: Transformation of each cell uses a weighted average of the two cells of 
# each anchor across anchors of the datasets. Weights determined by cell 
# similarity score (distance between cell and k nearest anchors) and anchor 
# scores, so cells in the same neighborhood should have similar correction values.

# If cell types are present in one dataset, but not the other, then the cells 
# will still appear as a separate sample-specific cluster.

# If sample has few cells integration will fail. So, we need to find kweight.
# Method1: Use half of minimum number of anchors between samples (NOT GOOD)
# anchors <- data.frame(dataset2=NA, n=NA)
# for (i in unique(integ_anchors@anchors$dataset1)){
#   test <- integ_anchors@anchors %>% 
#     dplyr::filter(dataset1==i) %>% 
#     dplyr::count(dataset2)
#   
#   anchors <- dplyr::bind_rows(anchors, test)
# }
# kweight1 <- floor(anchors %>%
#                     dplyr::select(n) %>%
#                     dplyr::arrange(n) %>% 
#                     dplyr::slice_head(n=1) %>% 
#                     unlist(., use.names = FALSE)/2)

#**STEP 8B: RUN PCA USING 3000 INTEGRATION FEATURES & UMAP USING FIRST 40 PCs**#
# You need to run PCA and UMAP after integration in order to visualize correctly
# because IntegrateData() uses a different set of 3000 variable genes. So, new
# PCs will need to be calculated.
# Note: If you used SCTransform() before integration, you don't need to run 
# ScaleData() after integration. However, if you ONLY used NormalizeData() 
# before integration, you need to use ScaleData() after integration.

integrate_data <- function(split_seurat, filtered_seurat){
  
  #*********STEP 7A: DECLARE REFERENCE SAMPLES FOR INTEGRATING THE DATA**********#
  # We assign index id, then select the index ids for appropriate reference samples
  if (is.na(ref1) & is.na(ref2)){
    ref_samples <- NULL
  } else if(!is.na(ref1) & is.na(ref2)){
    ref_samples <- as.vector(filtered_seurat@meta.data %>% 
                               dplyr::distinct_at("Sample", .keep_all = TRUE) %>%
                               dplyr::filter(!(Sample %in% samples_with_few_cells)) %>%
                               dplyr::mutate(index = row_number()) %>%
                               dplyr::filter(!!rlang::sym(ref1) == ref1_value) %>%
                               dplyr::select(index))[[1]]
  } else{
    ref_samples <- as.vector(filtered_seurat@meta.data %>% 
                               dplyr::distinct_at("Sample", .keep_all = TRUE) %>%
                               dplyr::filter(!(Sample %in% samples_with_few_cells)) %>%
                               dplyr::mutate(index = row_number()) %>%
                               dplyr::filter(!!rlang::sym(ref1) == ref1_value, !!rlang::sym(ref2) == ref2_value) %>%
                               dplyr::select(index))[[1]]
  }  
  
  cat("\nReference samples are:")
  for (i in ref_samples){
    cat("\n", i, ":", unique(filtered_seurat@meta.data$Sample)[i])  
  }
  
  #***STEP 7B: SELECT 3000 MOST VARIABLE GENES TO USE FOR INTEGRATING THE DATA***#
  integ_features <- Seurat::SelectIntegrationFeatures(object.list = split_seurat,
                                                      nfeatures = 3000,
                                                      assay = NULL, #c("SCT", "SCT"),
                                                      fvf.nfeatures = 3000)
  
  #******************STEP 7C: FIND RESIDUALS FOR MISSING GENES*******************#
  split_seurat <- Seurat::PrepSCTIntegration(object.list = split_seurat,
                                             assay = "SCT",
                                             anchor.features = integ_features,
                                             sct.clip.range = NULL)
  
  #******STEP 7D: FIND COMMON ANCHORS BETWEEN SAMPLES TO INTEGRATE THE DATA******#
  integ_anchors <- Seurat::FindIntegrationAnchors(object.list = split_seurat,
                                                  assay = NULL, #c("SCT", "SCT"),
                                                  reference = ref_samples,
                                                  anchor.features = integ_features,
                                                  scale = TRUE,
                                                  normalization.method = "SCT",
                                                  sct.clip.range = NULL,
                                                  reduction = "rpca",
                                                  l2.norm = TRUE,
                                                  dims = 1:30,
                                                  k.anchor = 5,
                                                  k.filter = 200,
                                                  k.score = 30,
                                                  max.features = 200,
                                                  nn.method = "annoy",
                                                  n.trees = 50,
                                                  eps = 0)
  
  #******STEP 7E: FIND OPTIMUM k.weight FOR USE IN Seurat::IntegrateData()*****#
  # Find minimum anchors between 2 datasets
  kweight1 <- as.data.frame(integ_anchors@anchors) %>%
    dplyr::group_by(dataset1, dataset2) %>%
    distinct_at("cell1", .keep_all=TRUE) %>%
    dplyr::summarize(n = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(n) %>% 
    unlist(use.names=FALSE) %>% 
    min()
 
  # Find half of number of cells in sample with least cell count
  kweight2 <- filtered_seurat@meta.data %>% 
    dplyr::count(Sample) %>% 
    dplyr::filter(n >=50) %>%
    dplyr::select(n) %>% 
    unlist(use.names=FALSE) %>% 
    min()
  
  kweight2 <- floor(kweight2/2)
    
  kweight <- base::min(kweight1, kweight2, 100)
  dplyr::if_else(kweight >= 100, 100, kweight)
  cat("\n", celltype, "\tkweight1:", kweight1, "\tkweight2:", kweight2, "\tkweight:",  kweight, "\n")
  
  # NOTE: Integration will not fail anymore. If it fails, identify the 2 
  # datasets that are involved in the error and use kweight = number of anchors 
  # for these 2 datasets.
  cat("\nNumber of unique anchors between datasets\n")
  print(as.data.frame(integ_anchors@anchors) %>% 
          dplyr::group_by(dataset1, dataset2) %>% 
          distinct_at("cell1", .keep_all=TRUE) %>% 
          dplyr::summarize(n = n()), n=1000)
  
  #************************STEP 8A: INTEGRATE THE DATA*************************#
  # NOTE: weight.reduction = NULL means new PCA will be calculated & used to
  # calculate anchor weights
  integrated_seurat <- Seurat::IntegrateData(anchorset = integ_anchors,
                                             new.assay.name = "integrated",
                                             normalization.method = "SCT",
                                             features = NULL,
                                             features.to.integrate = NULL,
                                             dims = 1:30,
                                             k.weight = kweight, #default is 100
                                             weight.reduction = NULL,
                                             sd.weight = 1)
  
  #**STEP 8B: RUN PCA USING 3000 INTEGRATION FEATURES & UMAP USING FIRST 40 PCs**#
  integrated_seurat <- Seurat::RunPCA(object = integrated_seurat,
                                      assay = "integrated",
                                      features = NULL)
  
  integrated_seurat <- Seurat::RunUMAP(object = integrated_seurat,
                                       dims = 1:40,
                                       reduction = "pca")
  
  return(integrated_seurat)
}

#******************************************************************************#
#                 CLUSTER THE CELLS & REMOVE SCARCE CLUSTERS                   #
#******************************************************************************#

# FindClusters performs graph-based clustering on the neighbor graph that is 
# constructed with the FindNeighbors function call. This neighbor graph is 
# constructed using PCA space when you specify reduction = "pca" in 
# FindNeighbors().

# NOTE: It is recommended to adjust k.param of FindNeighbors() [default=20] to 
# the same value as n.neighbors of UMAP() [default=30] 

cluster_data <- function(integrated_seurat){
  
  #****************STEP 9A: FIND NEAREST NEIGHBORS FOR EVERY CELL****************#
  # Determine the K-nearest neighbor graph
  integrated_seurat <- Seurat::FindNeighbors(object = integrated_seurat,
                                             reduction = "pca",
                                             dims = 1:40,
                                             assay = NULL,
                                             features = NULL,
                                             k.param = 20, #See NOTE
                                             return.neighbor = FALSE,
                                             compute.SNN = TRUE,
                                             prune.SNN = 1/15,
                                             nn.method = "annoy",
                                             n.trees = 50,
                                             annoy.metric = "euclidean",
                                             nn.eps = 0,
                                             force.recalc = FALSE,
                                             do.plot = FALSE,
                                             graph.name = NULL,
                                             l2.norm = FALSE,
                                             cache.index = FALSE)
  
  #*******STEP 9B: SEPARATE CELLS INTO CLUSTERS BASED ON NEAREST NEIGHBORS*******#
  # Determine the clusters for various resolutions
  integrated_seurat <- Seurat::FindClusters(object = integrated_seurat,
                                            resolution =c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                                            algorithm = 3,      #4 = Leiden
                                            method = "matrix",
                                            graph.name = NULL,
                                            modularity.fxn = 1,
                                            n.start = 10,
                                            n.iter = 10,
                                            random.seed = 0,
                                            group.singletons = TRUE)
  
  return(integrated_seurat)
}

# Remove clusters containing <= 3 cells so that UMAP is not littered with 80+ 
# cluster numbers.
clean_data <- function(integrated_seurat, active_res){
  
  # Identify clusters which have less than 3 cells
  sparse_clusters <- integrated_seurat@meta.data %>%
    dplyr::count(get(active_res)) %>%
    dplyr::filter(n <= 3) %>%
    dplyr::select(identity(1)) %>%
    unlist(.,use.names=FALSE) %>%
    as.character() %>%
    as.numeric()
  
  sparse_cells <- integrated_seurat@meta.data %>%
    dplyr::filter(get(active_res) %in% sparse_clusters) %>%
    dplyr::select(Cell) %>%
    unlist(.,use.names=FALSE)
  
  
  
  integrated_seurat <- subset(x = integrated_seurat,
                              Cell %in% sparse_cells,
                              invert = TRUE)
  
  cat("\nCells removed:", length(sparse_cells))
  
  return(integrated_seurat) 
}

save_data <- function(integrated_seurat){
  
  saveRDS(integrated_seurat, paste0(seurat_results, "integrated_seurat_snn", 
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
}

#******************************************************************************#
#                       VISUALIZE DATA AFTER INTEGRATION                       #
#******************************************************************************#

plot_post_integration <- function(integrated_seurat, active_res){
  
  col_id <- colnames(integrated_seurat@meta.data %>% 
                       dplyr::select(starts_with("integrated_snn_res")))
  
  # Create a function to plot PCA
  PCA_plot <- function(idents, split_by){
    
    # Set identity to an existing column in meta data
    Idents(object = integrated_seurat) <- idents
    
    if (is.na(split_by)){
      split_by <- NULL
    }
    
    Seurat::DimPlot(object = integrated_seurat,
                    reduction = "pca",
                    split.by = split_by,
                    label = TRUE,
                    raster = FALSE,
                    ncol = dplyr::if_else(length(unique(integrated_seurat@meta.data$Sample)) < 3, 2, 3),
                    combine = TRUE) +
      NoLegend() +
      ggplot2::labs(title = base::gsub(pattern="integrated_snn_res.", replacement="", x = idents))
    #facet_wrap(.~Sample, nrow = 4)
  }
  
  # Create a function to plot UMAP
  UMAP_plot <- function(idents, split_by){
    
    # Set identity to an existing column in meta data
    Idents(object = integrated_seurat) <- idents
    
    if (is.na(split_by)){
      split_by <- NULL
    }
    
    Seurat::DimPlot(object = integrated_seurat,
                    reduction = "umap",
                    split.by = split_by,
                    label = TRUE,
                    raster = FALSE,
                    ncol = dplyr::if_else(length(unique(integrated_seurat@meta.data$Sample)) < 3, 2, 3),
                    combine = TRUE) +
      NoLegend() +
      ggplot2::labs(title = base::gsub(pattern="integrated_snn_res.", replacement="", x = idents))
    #facet_wrap(.~Sample, nrow = 4)  
  }
  
  # Create a function to plot numerical metrics such as nUMIs/cell, nGenes/cell, 
  # S phase, G2M phase, CC scores and mitoRatio using UMAP
  UMAP_plot_numerical_metrics <- function(idents){
    
    # Set identity to an existing column in meta data
    Idents(object = integrated_seurat) <- idents
    
    Seurat::FeaturePlot(object = integrated_seurat,
                        features = c("nUMIs", "nGenes", "S.Score", "G2M.Score", 
                                     "CC.Score", "MitoRatio"),
                        min.cutoff = 'q10',
                        pt.size = 0.4,
                        reduction = "umap",
                        order = TRUE,  # plot positive cells above negative cells
                        label = TRUE,
                        raster = FALSE,
                        ncol = 2,
                        combine = TRUE)
  }
  
  # Create a function to plot categorical metrics such as Sample, Phase, 
  # Condition using DimPlot() as these cannot be plotted using FeaturePlot()
  UMAP_plot_categorical_metrics <- function(idents, split_by){
    
    # Set identity to an existing column in meta data
    Idents(object = integrated_seurat) <- idents
    
    if (is.na(split_by)){
      split_by <- NULL
    }
    
    Seurat::DimPlot(object = integrated_seurat,
                    split.by = split_by,
                    pt.size = 0.4,
                    reduction = "umap",
                    order = TRUE,  # plot positive cells above negative cells
                    label = FALSE,
                    raster = FALSE,
                    ncol = dplyr::if_else(split_by == "Sample", 4, 
                                          dplyr::if_else(split_by == "Phase", 3, 2)),
                    combine = TRUE) +
      NoLegend()    #facet_wrap(.~metrics[i], nrow = 4)
  }
  
  #****************STEP 10A: PLOT PCA & UMAP AFTER INTEGRATION*****************#
  #******STEP 10B: PLOT NUMERICAL & CATEGORICAL METRICS AFTER INTEGRATION******#
  #******************STEP 10C: PLOT UMAP FOR EVERY RESOLUTION******************#
  
  funcs <- c("PCA_plot", "UMAP_plot", "UMAP_plot_numerical_metrics", 
             "UMAP_plot_categorical_metrics", "UMAP_plot")
  
  for (j in 1:length(funcs)){
    
    if (j == 1){
      filename <- paste0("PCA_plot_Post_Integration_", celltype, ".jpg")
      get(funcs[j])(active_res, "Sample")
    } else if (j == 2){
      filename <- paste0("UMAP_plot_Post_Integration_", celltype, ".jpg")
      get(funcs[j])(active_res, "Sample")
    } else if (j == 3){
      filename <- paste0("UMAP_plot_numerical_metrics_", celltype, ".jpg")
      get(funcs[j])(active_res)
    } else if (j == 4){
      filename <- paste0("UMAP_plot_categorical_metrics_", celltype, ".jpg")
      purrr::map2(.x = c(active_res, active_res, active_res),
                  .y = c("Sample", "Phase", "Condition"),
                  .f = UMAP_plot_categorical_metrics) %>%
        cowplot::plot_grid(plotlist = .,
                           align = "hv",
                           axis = "tblr",
                           nrow = NULL,
                           ncol = 3,
                           rel_widths = 1,
                           rel_heights = 1,
                           greedy = TRUE,
                           byrow = TRUE)
    } else if (j == 5){
      filename <- paste0("UMAP_all_resolutions_", celltype, ".jpg")
      purrr::map2(.x = col_id,
                  .y = base::rep(x = NA, times = length(col_id)),
                  .f = UMAP_plot) %>%
        cowplot::plot_grid(plotlist = .,
                           align = "hv",
                           axis = "tblr",
                           nrow = NULL,
                           ncol = 3,
                           rel_widths = 1,
                           rel_heights = 1,
                           greedy = TRUE,
                           byrow = TRUE)
    }
    
    # Save the plot
    ggplot2::ggsave(filename = filename,
                    plot = last_plot(),
                    device = "jpeg",
                    path = diagnostics_path,
                    scale = 1,
                    width = dplyr::case_when(j == 1 | j==2 | j==3 ~ 8.5,
                                             j == 4 ~ 8.5*3,
                                             j ==5 ~ 8.5*2),
                    height = 11,
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = "white")
  }
}

#******************************************************************************#
#                            CLUSTER IDENTIFICATION                            #
#******************************************************************************#

# NOTE: Perform cluster identification using conserved markers. If you cannot
# identify any cluster using conserved markers, then use the markers from
# FindAllMarkers() to identify ONLY the unidentified clusters.
# FindConservedMarkers() is more accurate as it gives genes conserved in
# multiple samples. FindAllMarkers() merges all cluster 0 cells from all samples
# and compares it to rest of clusters (all the cells from all samples)

plot_modules <- function(active_res, sheetname){
  
  #**********STEP 11A: PLOT MODULE SCORES ON UMAP TO IDENTIFY CELL TYPES*********#
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Set "RNA" assay to be the default assay.
  DefaultAssay(integrated_seurat) <- "RNA"
  
  # Change active.ident to major clusters we annotated
  Idents(object = integrated_seurat) <- active_res
  
  # Read top markers
  top_markers <- openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, "KMS_Markers.xlsx"),
                                     sheet = sheetname)
  
  # Plot modules of top markers on UMAP to identify clusters
  module_plot <- function(i){
    
    feature_subset <- as.list(top_markers %>%
                                dplyr::filter(CELL_TYPE == i) %>%
                                dplyr::select(dplyr::if_else(species == "Homo sapiens", "HUMAN_GENE", "MOUSE_GENE")))[[1]]
    
    feature_subset <- intersect(feature_subset, rownames(integrated_seurat@assays$RNA@counts))
    
    if (length(feature_subset) > 0){
      x <- Seurat::AddModuleScore(object = integrated_seurat,
                                  features = list(feature_subset),
                                  name = make.names(i),
                                  pool = NULL,
                                  nbin = 24,
                                  ctrl = 5, #100,
                                  k = FALSE,
                                  assay = "RNA",
                                  seed = 1,
                                  search = FALSE)
      
      Seurat::FeaturePlot(object = x,
                          slot = "data",
                          features = paste0(make.names(i),1),
                          cols =  c("grey", viridis(n = 10, option = "C", direction = -1)),
                          pt.size = 0.4,
                          order = TRUE,
                          min.cutoff = 'q10',
                          reduction = "umap",
                          label = TRUE,
                          combine = TRUE,
                          raster = FALSE)
    }
  }
  
  c <- unique(top_markers$CELL_TYPE)
  # Plot UMAPs
  purrr::map(.x = c, 
             .f = module_plot) %>%
    cowplot::plot_grid(plotlist = .,
                       align = "hv",
                       axis = "tblr",
                       nrow = NULL,
                       ncol = dplyr::if_else(length(c) > 10, 5, ceiling(length(c)/3)),
                       rel_widths = 1,
                       rel_heights = 1,
                       greedy = TRUE,
                       byrow = TRUE)
  
  ggplot2::ggsave(filename = paste0("Module_plot(", sheetname, ")_", celltype, ".jpg"),
                  plot = last_plot(),
                  device = "jpeg",
                  path = diagnostics_path,
                  scale = 1,
                  width = 8.5*4,
                  height = 11*2,
                  units = c("in"),
                  dpi = 300,
                  limitsize = FALSE,
                  bg = "white")
}

plot_genes <- function(active_res, sheetname){
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Set "RNA" assay to be the default assay.
  DefaultAssay(integrated_seurat) <- "RNA"
  
  # Change active.ident to major clusters we annotated
  Idents(object = integrated_seurat) <- active_res
  
  # Read markers
  features <- openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, "KMS_Markers.xlsx"),
                                  sheet = sheetname)
  
  features <- as.list(features %>%
                        dplyr::select(dplyr::if_else(species == "Homo sapiens", "HUMAN_GENE", "MOUSE_GENE")))[[1]]
  
  features <- intersect(features, rownames(integrated_seurat@assays$RNA@data))
  
  feature_plot <- function(i){
    
    split <- "Condition" #NULL
    Seurat::FeaturePlot(object = integrated_seurat,
                        slot = "data",
                        features = i,
                        split.by = split,
                        cols =  c("grey", viridis(n = 10, option = "C", direction = -1)),
                        pt.size = 0.4,
                        order = TRUE,
                        min.cutoff = 'q10',
                        reduction = "umap",
                        label = TRUE,
                        combine = TRUE,
                        raster = FALSE)
  }
  
  # Plot UMAPs
  for (n in 1:ceiling((length(features)/12))){
    j <- 12*n-11
    k <- 12*n
    vec <- features[j:k]
    purrr::map(.x = vec[!is.na(vec)], 
               .f = feature_plot) %>%
      cowplot::plot_grid(plotlist = .,
                         align = "hv",
                         axis = "tblr",
                         nrow = 2,
                         ncol = 6,
                         rel_widths = 1,
                         rel_heights = 1,
                         greedy = TRUE,
                         byrow = TRUE)
    
    ggplot2::ggsave(filename = paste0("Feature_plot_",n, ".jpg"),
                    plot = last_plot(),
                    device = "jpeg",
                    #path = diagnostics_path,
                    scale = 1,
                    width = dplyr::if_else(is.null(split), 9*3, 9*6),
                    height = dplyr::if_else(is.null(split), 11, 11),
                    units = c("in"),
                    dpi = 300,
                    limitsize = FALSE,
                    bg = "white")
  }
}

remove_mixed <- function(active_res, Mixed){
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  Idents(integrated_seurat) <- active_res
  cat("\n", celltype, " present initially:", nrow(integrated_seurat@meta.data))
  
  integrated_seurat <- subset(x = integrated_seurat,
                              !!rlang::sym(active_res) %in% Mixed[[celltype]],
                              invert = TRUE)
  cat("\n", celltype, " present finally:", nrow(integrated_seurat@meta.data))
  
  saveRDS(integrated_seurat, paste0(seurat_results, "integrated_seurat_snn", 
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
}

get_markers <- function(annotations, active_res){
  
  #****************STEP 11B: IDENTIFY ALL MARKERS FOR EACH CLUSTER***************#
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Set "RNA" assay to be the default assay.
  DefaultAssay(integrated_seurat) <- "RNA"
  
  # Change active.ident to major clusters we annotated
  Idents(object = integrated_seurat) <- active_res
  
  # Get all markers for any given cluster
  all_markers <- Seurat::FindAllMarkers(object=integrated_seurat,
                                        assay = "RNA",
                                        features = NULL,
                                        logfc.threshold = 0.25,
                                        test.use = "wilcox",
                                        slot = "data",
                                        min.pct = 0.1,
                                        min.diff.pct = 0.1,
                                        node = NULL,
                                        verbose = TRUE,
                                        only.pos = TRUE,
                                        max.cells.per.ident = Inf,
                                        random.seed = 1,
                                        latent.vars = NULL,
                                        min.cells.feature = 3,
                                        min.cells.group = 1,
                                        pseudocount.use = 1,
                                        mean.fxn = NULL,
                                        fc.name = NULL,
                                        base = 2,
                                        return.thresh = 0.01,
                                        densify = FALSE)
  
  all_markers <- all_markers %>% 
    dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio = pct.1/pct.2) %>%
    dplyr::left_join(y = unique(annotations[, c("SYMBOL", "CHR", "DESCRIPTION")]), by = c("gene" = "SYMBOL")) %>%
    dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION)
  
  # Find top markers for each major cluster
  top_markers <- all_markers %>%
    dplyr::filter(avg_log2FC >= 0.58 & p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(ratio), desc(avg_log2FC)) %>%
    dplyr::slice_head(n=20) %>%
    ungroup()
  
  # Save all the markers
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "All_Markers")
  openxlsx::writeData(wb = wb, sheet = "All_Markers", x = all_markers)
  openxlsx::addWorksheet(wb = wb, sheetName = "Top_Markers")
  openxlsx::writeData(wb = wb, sheet = "Top_Markers", x = top_markers)
  openxlsx::saveWorkbook(wb = wb, file = paste0(seurat_results, "Markers_All_", active_res, "_", celltype, ".xlsx"), overwrite = TRUE)
}  

#******************************************************************************#
#             STEP 11: SUBSET SPECIFIC CELLTYPE FOR SUBTYPE ANALYSIS           #
#******************************************************************************#

prep_data <- function(integrated_seurat, celltype){
  
  cells <- unique(integrated_seurat@meta.data$cell_type)
  cells <- cells[grepl(pattern = celltype, x = cells)]
  
  # Keep ONLY necessary celltype being analyzed
  filtered_seurat <- subset(x = integrated_seurat, 
                            cell_type %in% cells)
  
  # Print cell numbers to double check
  print(filtered_seurat@meta.data %>% dplyr::count(Condition, cell_type, sub_type, Sample))
  
  # Remove unwanted assays after changing default assay
  DefaultAssay(filtered_seurat) <- "RNA"
  filtered_seurat[["SCT"]] <- NULL
  filtered_seurat[["integrated"]] <- NULL
  
  return(filtered_seurat)
}

#******************************************************************************#
#                      STEP 12: IDENTIFY SIMILAR CLUSTERS                      # 
#******************************************************************************#

identify_lineage <- function(celltype){
  
  cat("\n**********************", celltype, "Analysis**********************\n")
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Create matrix similar to read data for each cluster to identify closely related clusters
  cor_meta_data <- integrated_seurat@meta.data %>%
    dplyr::distinct(integrated_snn_res.1.4, .keep_all = TRUE) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("integrated_snn_res.1.4")
  
  clusters <- as.vector(unlist(integrated_seurat@meta.data %>% dplyr::distinct(integrated_snn_res.1.4)))
  cor_read_data <- matrix(NA, nrow = nrow(integrated_seurat@assays$RNA@counts), ncol = length(clusters))
  rownames(cor_read_data) <- rownames(integrated_seurat@assays$RNA@counts)
  colnames(cor_read_data) <- clusters
  
  for(i in 1:length(clusters)){
    
    # Create a list of cells for each cluster
    cor_cells_subset <- rownames(integrated_seurat@meta.data %>% dplyr::filter(integrated_snn_res.1.4 == as.numeric(clusters[i])))
    
    # Get the raw counts of these cells. counts are stored in sparse matrix format
    # So, use data.frame to convert "." in sparse matrix to "0"
    cor_subset <- data.frame(integrated_seurat@assays$RNA@counts[,cor_cells_subset])
    cor_read_data[,i]  <- rowSums(cor_subset)
  }
  
  # Remove unwanted genes
  unwanted_genes <- rownames(cor_read_data)[grep(pattern = "^RP[SL]|RIK$|^MT-|^GM[0-9.]+$", x= rownames(cor_read_data), ignore.case=TRUE)]
  keep_genes <- setdiff(rownames(cor_read_data), unwanted_genes)  
  cor_read_data <- cor_read_data[keep_genes,]
  
  # Create DESeq2 object so we can use DESeq2::rlog() or DESeq2::vst() on it
  cor_dds <- DESeqDataSetFromMatrix(countData = cor_read_data,
                                    colData = cor_meta_data,
                                    design = ~1)
  
  #Determine number of cuts in heatmap
  cuts <- ceiling(length(clusters)/3)
  
  # Use vst as well as rld
  vst <- DESeq2::vst(cor_dds)
  head(assay(vst))
  sampleDists <- dist(t(assay(vst)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  pheatmap::pheatmap(mat = sampleDistMatrix,
                     scale = "none",
                     cutree_rows = cuts,
                     cutree_cols = cuts,
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = TRUE,   #cluster the columns
                     fontsize = 8, 
                     fontsize_row = 8, 
                     fontsize_col = 8,
                     angle_col = c("270", "0", "45", "90", "315"),
                     fontsize_number = 0.8*fontsize,
                     filename = paste0(diagnostics_path, "Cluster_Correlation_vst_", celltype, ".jpg"))
  
  rlog <- DESeq2::rlog(cor_dds)
  head(assay(rlog))
  sampleDists <- dist(t(assay(rlog)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  pheatmap::pheatmap(mat = sampleDistMatrix,
                     scale = "none",
                     cutree_rows = cuts,
                     cutree_cols = cuts,
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = TRUE,   #cluster the columns
                     fontsize = 8, 
                     fontsize_row = 8, 
                     fontsize_col = 8,
                     angle_col = c("270", "0", "45", "90", "315"),
                     fontsize_number = 0.8*fontsize,
                     filename = paste0(diagnostics_path, "Cluster_Correlation_rlog_", celltype, ".jpg"))
}


#******************************************************************************#
#       RE-ANNOTATE ORIGINAL SEURAT OBJECT WITH SUBTYPE, CELLCLASS INFO        #
#******************************************************************************#

re_annotate <- function(celltypes){
  
  celltype <- celltypes[1]
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  merged_metadata <- integrated_seurat@meta.data
  
  
  for (celltype in celltypes[-1]){
    
    # Load the integrated seurat object
    integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                        dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
    
    merged_metadata <- dplyr::bind_rows(merged_metadata, integrated_seurat@meta.data)
  }
  
  # Load the full integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))
  
  metadata <- integrated_seurat@meta.data

  # You can see there are more cells in metadata than merged_metadata.   
  # Some of these cells were likely removed during subtype analysis while others
  # cells belong to cell_type like Muscle, Neurons, Unclassified. So, we need to
  # generate sub_type and cell_class info for these cells.
  missing_metadata <- dplyr::anti_join(x = metadata,
                                       y = merged_metadata, 
                                       by=c("Cell"="Cell")) %>%
    dplyr::mutate(cell_class = dplyr::case_when(cell_type %in% base::setdiff(x=unique(metadata$cell_type), y=unique(merged_metadata$cell_type)) ~ cell_type,
                                                cell_type == "Unclassified" ~ "Unclassified",
                                                TRUE ~ "Mixed"),
                  sub_type = dplyr::case_when(cell_type %in% base::setdiff(x=unique(metadata$cell_type), y=unique(merged_metadata$cell_type)) ~ cell_type,
                                              cell_type == "Unclassified" ~ "Unclassified",
                                              TRUE ~ "Mixed"))
  
  # Import cell_class, sub_type and cell_type from merged_metadata
  subtyped_metadata <-  dplyr::inner_join(x = metadata %>% dplyr::select(everything(), -cell_type, -sub_type, -cell_class),
                                          y = merged_metadata %>% dplyr::select(Cell, cell_type, sub_type, cell_class), 
                                          by=c("Cell"="Cell"))
  
  final_metadata <- dplyr::bind_rows(subtyped_metadata, missing_metadata)
  rownames(final_metadata) <- final_metadata$Cell   #this is most important step, else UMAP labels will be wrong
  
  # Print counts before annotation
  integrated_seurat@meta.data %>% count(cell_type, sub_type, cell_class)
  
  # Re-annotate
  integrated_seurat@meta.data <- final_metadata
  
  # Print counts after annotation
  integrated_seurat@meta.data %>% count(cell_type,sub_type, cell_class)
  
  # # Remove cross labelled cells and cell_class to metadata
  # # Some cells may have cell_type "Myeloid - MDSC" but sub_type as "Myeloid- cDcs".
  # integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
  #   dplyr::mutate(cell_class = dplyr::if_else(cell_type == sub_type |
  #                                               cell_type == "Lymphoid - B" & grepl(pattern = "B|Plasma", x=sub_type) |
  #                                               cell_type == "Lymphoid - T" & grepl(pattern = "CD4|CD8|Gamma|NKT", x=sub_type) |
  #                                               cell_type == "Myeloid - Macrophages, DCs" & grepl(pattern = "Macrophage|cDCs|pDCs", x=sub_type) |
  #                                               cell_type == "Epithelial" & grepl(pattern = "Epithelial", x=sub_type) |
  #                                               cell_type == "Fibroblasts" & grepl(pattern = "Fibroblasts", x=sub_type),
  #                                             gsub(pattern = "\ -.*",replacement = "",x = cell_type), "Mixed"))
  
  # Save the object after checking UMAP
  saveRDS(integrated_seurat, paste0(seurat_results, "integrated_seurat_snn.rds"))
}

#******************************************************************************#
#              VISUALIZE CELLTYPES, SUBTYPES, CELLCLASSES IN UMAP              #
#******************************************************************************#

# You must define celltype and split as global variables
visualize_UMAP <- function(){
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Remove unwanted cells and samples
  plot_object <- subset(x = integrated_seurat,
                        subset = (cell_class %in% c("Mixed", "Unclassified")),
                        invert = TRUE)
  
  gg1 <- Seurat::DimPlot(object = integrated_seurat,
                         reduction = "umap",
                         cols = my_palette,
                         label = FALSE,
                         group.by = "cell_class",
                         split.by = split,
                         shape.by = NULL,
                         pt.size = 0.2,
                         label.size = 5,
                         repel = FALSE,
                         raster = FALSE) +
    ggplot2::labs(fill = "CLUSTERS",
                  x = "UMAP_1",
                  y = "UMAP_2") +
    my_theme
  
  gg2 <- Seurat::DimPlot(object = plot_object,
                         reduction = "umap",
                         cols = my_palette,
                         label = FALSE,
                         group.by = "cell_type",
                         split.by = split,
                         shape.by = NULL,
                         pt.size = 0.2,
                         label.size = 5,
                         repel = FALSE,
                         raster = FALSE) +
    ggplot2::labs(fill = "CLUSTERS",
                  x = "UMAP_1",
                  y = "UMAP_2") +
    my_theme
  
  gg3 <- Seurat::DimPlot(object = plot_object,
                         reduction = "umap",
                         cols = my_palette,
                         label = FALSE,
                         group.by = "sub_type",
                         split.by = split,
                         shape.by = NULL,
                         pt.size = 0.2,
                         label.size = 5,
                         repel = FALSE,
                         raster = FALSE) +
    ggplot2::labs(fill = "CLUSTERS",
                  x = "UMAP_1",
                  y = "UMAP_2") +
    my_theme
  
  ggplot2::ggsave(filename = paste0("UMAP_", celltype, "_final_", split, ".pdf"),
                  plot = gg1+gg2+gg3,
                  device = "pdf",
                  path = seurat_results,
                  scale = 1,
                  width = dplyr::if_else(is.null(split), 9*3, 9*2),
                  height = dplyr::if_else(is.null(split), 11, 11*3),
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                  VERIFY CLUSTER ANNOTATION USING DOT PLOTS                   #
#******************************************************************************#

# Excel file MUST contain column (i) CELL_TYPE (ii) HUMAN_GENE (iii) MOUSE_GENE
# Excel file MUST be named KMS_Markers.xlsx and located in scripts_path
# Excel file MUST have sheetnames identical to celltype variable (global)
visualize_dotplot <- function(){
  
  # Read file containing marker genes
  kms_markers <- openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, "KMS_Markers.xlsx"),
                                     sheet = celltype)
  
  # Remove all duplicated genes
  feature <- kms_markers %>%
    dplyr::select(dplyr::if_else(species == "Homo sapiens", "HUMAN_GENE", "MOUSE_GENE")) %>%
    dplyr::distinct() %>%
    unlist(use.names = FALSE)
  
  feature <- feature[!is.na(feature)]
  
  if (celltype == "All Markers"){
    celltype <- NULL
  }
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Remove unwanted cells and samples
  integrated_seurat <- subset(x = integrated_seurat,
                              subset = (cell_class %in% c("Mixed", "Unclassified")),
                              invert = TRUE)
  
  # Set Idents to cell_type so UMAPs show the Idents instead of cluster number
  if (is.null(celltype)){
    Idents(integrated_seurat) <- "cell_type"
  } else {
    Idents(integrated_seurat) <- "sub_type"
  }
  
  # We re-order the active ident alphabetically
  Idents(integrated_seurat) <- base::factor(x = integrated_seurat@active.ident, 
                                            levels = sort(levels(integrated_seurat@active.ident)))
  
  plot_dot <- function(i){
    feature_subset <- feature[((i*20)-19):(i*20)]
    feature_subset <- feature_subset[!is.na(feature_subset)]
    
    Seurat::DotPlot(object = integrated_seurat,
                    assay = "RNA",
                    features = feature_subset,
                    #cols = c("blue", "red"),
                    #col.min = -2,
                    #col.max = 2,
                    dot.min = 0,
                    dot.scale = 4,
                    idents = NULL,
                    group.by = NULL,
                    split.by = NULL,
                    cluster.idents = FALSE,
                    scale = TRUE,
                    scale.by = "size",
                    scale.min = 0,
                    scale.max = 100) +
      ggplot2::geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.25) + #stroke is width of circle
      ggplot2::scale_colour_distiller(type = "div", palette = "RdYlGn", direction = -1) +
      ggplot2::guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white", stroke=0.75))) +
      my_theme
  }
  
  purrr::map(.x = 1:ceiling(length(feature)/20), 
             .f = plot_dot) %>%
    cowplot::plot_grid(plotlist = .,
                       align = "hv",
                       axis = "tblr",
                       nrow = ceiling(length(feature)/20),
                       ncol = 1,
                       rel_widths = 1,
                       rel_heights = 1,
                       greedy = TRUE,
                       byrow = TRUE)
  
  ggplot2::ggsave(filename = paste0("Cluster_Verification_Dot_plot_", celltype, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = seurat_results,
                  scale = 1,
                  width = 11, #5+0.5*length(feature),
                  height = 4*ceiling(length(feature)/20), #1+0.5*length(levels(Idents(integrated_seurat))),
                  units = c("in"),
                  dpi = 300,
                  limitsize = FALSE,
                  bg = "white")
}

#******************************************************************************#
#                  STEP 13: CALCULATE STATS FOR EACH CLUSTER                   #               
#******************************************************************************#

#**STEP 13: CALCULATE NUMBER OF CELLS IN EACH CLUSTER FOR EVERY RESOLUTION**#
# NOTE: We identify sparse as well as dominant cells at this step but we 
# remove them after checking all plots in prior steps (10A-C)
calc_stats <- function(){
  wb <- openxlsx::createWorkbook()
  dominant_cells <- data.frame(Cell = NA)
  isolated_cells <- data.frame(Cell = NA)
  
  for (x in 1:length(col_id)){  
    
    # Calculate cells per cluster
    cells_per_cluster <- integrated_seurat@meta.data %>% 
      dplyr::group_by(!!rlang::sym(col_id[x]), Sample) %>%
      dplyr::count() %>%
      dplyr::rename(cluster = col_id[x], nCells = n) %>%
      tidyr::pivot_wider(id_cols = Sample, names_from = cluster, values_from = nCells) %>%
      base::replace(is.na(.), 0)
    
    # Calculate percent of cells in each cluster for each sample
    # Divide by rowSums to adjust for difference in cell number between samples
    cells_per_cluster_percent <- cells_per_cluster %>%
      dplyr::mutate(across(.cols = everything())*100/rowSums(across()))
    
    # Calculate total cells in each cluster across all samples
    cells_per_cluster_total <- c(list(Sample = "Total Cells per cluster"), colSums(cells_per_cluster[-1]))
    
    # Merge all data
    cells_per_cluster <- dplyr::bind_rows(cells_per_cluster, data.frame(data = " "), 
                                          cells_per_cluster_percent, data.frame(data = " "),
                                          cells_per_cluster_total) %>%
      dplyr::select(-data)
    
    # Save data
    openxlsx::addWorksheet(wb = wb, sheetName = paste0("Resolution_",all_res[x]))
    openxlsx::writeData(wb = wb, sheet = paste0("Resolution_",all_res[x]), x = cells_per_cluster)
    
    #****************************************************************************#
    # Set a value for outlier classification (outlier > 75% + lenient*IQR)
    # High lenient value means higher upper cutoff. So, filtering is lenient
    lenient <- 25
    
    # Calculate quartiles for each cluster excluding Sample column
    # Column names of transposed quartile  are:  0%, 25%,  50%,  75% and  100%
    # If value > 75% + 100*IQR, then we denote these clusters as outlier 
    quartiles <- apply(X = cells_per_cluster_percent[,-1], MARGIN = 2, FUN = stats::quantile) %>%
      t() %>%
      data.frame() %>%
      dplyr::mutate(upper_cutoff = .[[4]] + lenient*(.[[4]] - .[[2]])) %>%
      dplyr::select(upper_cutoff) %>%
      t()
    rownames(quartiles) <- NULL
    
    # Identify dominant cells
    for (i in 1:nrow(cells_per_cluster_percent)) {
      for (j in 1:ncol(quartiles)) {
        if (cells_per_cluster_percent[i, (j+1)] > quartiles[1, j]) {
          
          cluster_id <- colnames(cells_per_cluster_percent)[j+1]
          sample_id <- cells_per_cluster_percent$Sample[i]
          
          cells <- integrated_seurat@meta.data %>%
            dplyr::filter(!!rlang::sym(col_id[x]) == cluster_id & Sample == sample_id) %>%
            dplyr::select(Cell)
          
          dominant_cells <- rbind(dominant_cells, cells)
          
          cat("\nSample:", sample_id, 
              "\tCluster:", cluster_id,
              "\tCluster:", colnames(quartiles)[j],
              "\tcutoff:", quartiles[1,j])
        }
      }
    }
    
    # Identify sparse cells
    cells <- integrated_seurat@meta.data %>% 
      dplyr::add_count(get(col_id[x])) %>%
      dplyr::filter(n < 5) %>%
      dplyr::select(Cell)
    
    isolated_cells <- rbind(isolated_cells, cells)
    
    print(dim(isolated_cells))
  }
  
  # Identify cells consistently dominant in 3 of 6 resolutions)
  dominant_cells <- dominant_cells %>% 
    dplyr::count(Cell) %>%
    dplyr::filter(n >= 3)
  
  # Identify cells consistently sparsely clustered in all 6 resolutions
  isolated_cells <- isolated_cells %>% 
    dplyr::count(Cell) %>% 
    dplyr::filter(n == 6)
  
  # Merge
  if (nrow(dominant_cells) > 0 & nrow(isolated_cells) > 0){
    bad_cells <- dplyr::bind_rows(dominant_cells, isolated_cells) %>%
      dplyr::distinct_at("Cell")
  } else if (nrow(dominant_cells) > 0 & nrow(isolated_cells) == 0){
    bad_cells <- dominant_cells
  } else{
    bad_cells <- data.frame(Cell = "NA") 
  }
  
  # Save data
  openxlsx::addWorksheet(wb = wb, sheetName = paste0("Bad cells"))
  openxlsx::writeData(wb = wb, sheet = paste0("Bad cells"), x = bad_cells)
  openxlsx::saveWorkbook(wb = wb, file = paste0(seurat_results, "Cells_per_cluster_", celltype, ".xlsx"), overwrite = TRUE)
  
  # # Next, explore how well our clusters separate by the different PCs.
  # # To visualize this information, we need to extract the UMAP coordinate 
  # # information for the cells along with their corresponding scores for each of 
  # # the PCs to view by UMAP
  # 
  # for(i in 1:5){
  #   j <- 9*i-8
  #   k <- 9*i
  #   Seurat::FeaturePlot(object = integrated_seurat,
  #                       features = c(paste0("PC_",j:k)),
  #                       pt.size = 0.4,
  #                       order = TRUE,
  #                       min.cutoff = 'q10',
  #                       reduction = "umap",
  #                       label = TRUE,
  #                       combine = TRUE,
  #                       raster = FALSE)
  #   
  #   ggplot2::ggsave(filename = paste0("UMAP_for_PC_", j, "_through_", k, ".jpg"),
  #                   plot = last_plot(),
  #                   device = "jpeg",
  #                   path = diagnostics_path,
  #                   scale = 1,
  #                   width = 8.5,
  #                   height = 11,
  #                   units = c("in"),
  #                   dpi = 600,
  #                   limitsize = TRUE,
  #                   bg = "white")
  # }
  
  #*******STEP 10E: CHECK IF CELL CYCLE IS A SOURCE OF VARIATION VIA PCs*******#
  
  # # Save the top 30 (positive) and bottom 30 (negative) genes for each of 50 PCs
  # # Rows will be genes and columns will be PCs
  # pc_genes <- data.frame()[1:30, ]
  # rownames(pc_genes) <- paste0("Gene_",1:30)
  # for (i in 1:50){
  #   pc_genes[,2*i-1] <- names(sort(integrated_seurat@reductions$pca@feature.loadings[,i], decreasing = FALSE)[1:30])
  #   colnames(pc_genes)[2*i-1] <- paste0("Negative_PC_",i)
  #   pc_genes[,2*i] <- names(sort(integrated_seurat@reductions$pca@feature.loadings[,i], decreasing = TRUE)[1:30])
  #   colnames(pc_genes)[2*i] <- paste0("Positive_PC_",i)
  # }
  # 
  # # Find which PCs are affected by cell cycle genes
  # PCs_affected <- pc_genes
  # for(i in 1:ncol(PCs_affected)){
  #   PCs_affected[,i] <- pc_genes[,i] %in% cell_cycle_genes
  # }
  # 
  # # Populate the TRUE values with genes names
  # for(i in 1:ncol(PCs_affected)){
  #   for(j in 1:nrow(PCs_affected)){
  #     if(PCs_affected[j,i] == TRUE){
  #       PCs_affected[j,i] <- pc_genes[j,i]
  #     } else{
  #       PCs_affected[j,i] <- ""
  #     }
  #   }
  # }
  # 
  # # Remove PCs with no genes
  # PCs_affected <- PCs_affected[, colSums(PCs_affected != "") != 0]
  # 
  # # Save the dataframe in xlsx format
  # wb <- openxlsx::createWorkbook()
  # openxlsx::addWorksheet(wb = wb, sheetName = "PCA_genes")
  # openxlsx::writeData(wb = wb, sheet = "PCA_genes", x = pc_genes)
  # openxlsx::addWorksheet(wb = wb, sheetName = "Cell_cycle_affected_PCs")
  # openxlsx::writeData(wb = wb, sheet = "Cell_cycle_affected_PCs", x = PCs_affected)
  # openxlsx::saveWorkbook(wb = wb, file = paste0(diagnostics_path, "PCA_Genes_After_Integration.xlsx"), overwrite = TRUE)
  
  #*****STEP 10F: CHECK IF CELL CYCLE IS A SOURCE OF VARIATION VIA GRAPHS******#
  
  # integrated_seurat_cc <- RunPCA(object = integrated_seurat,
  #                                assay = NULL,
  #                                features = c(s_genes, g2m_genes),
  #                                npcs = 50,
  #                                rev.pca = FALSE,
  #                                weight.by.var = TRUE,
  #                                verbose = TRUE,
  #                                ndims.print = 1:5,
  #                                nfeatures.print = 30,
  #                                reduction.name = "pca",
  #                                reduction.key = "PC_",
  #                                seed.use = 42)
  # 
  # # Plot the PCA colored by cell cycle phase
  # Seurat::DimPlot(object = integrated_seurat_cc, 
  #                 reduction = "pca", 
  #                 group.by = "Phase",
  #                 split.by = "Sample",
  #                 raster = FALSE,
  #                 ncol = dplyr::if_else(length(unique(integrated_seurat_cc@meta.data$Sample)) < 3, 2, 3),
  #                 combine = TRUE) +
  #   #facet_wrap(.~Sample, nrow = 4) +				
  #   NoLegend()  
  # 
  # # Save the plot
  # ggplot2::ggsave(filename = "PCA_using_cell_cycle_genes.jpg",
  #                 plot = last_plot(),
  #                 device = "jpeg",
  #                 path = diagnostics_path,
  #                 scale = 1,
  #                 width = 8.5,
  #                 height = 11,
  #                 units = c("in"),
  #                 dpi = 600,
  #                 limitsize = TRUE,
  #                 bg = "white")
  
  #********STEP 10G: DETERMINE THE 'DIMENSIONALITY' FROM HEATMAP OF PCS********#
  
  # # cells parameter in DimHeatmap() specifies the number of cells with most 
  # # negative or postive PCA scores to be used for plotting. The idea is that we 
  # # are looking for a PC where the heatmap starts to look more “fuzzy”, i.e. where
  # # the distinctions between the groups of genes is not so distinct.
  # for(i in 1:5){
  #   j <- 9*i-8
  #   k <- 9*i
  #   Seurat::DimHeatmap(object = integrated_seurat,
  #                      dims = j:k,
  #                      nfeatures = 30,
  #                      cells = 1000,
  #                      reduction = "pca",
  #                      disp.min = -2.5,
  #                      disp.max = NULL,
  #                      balanced = TRUE,
  #                      projected = FALSE,
  #                      ncol = 3,
  #                      fast = FALSE,
  #                      raster = FALSE,
  #                      slot = "scale.data",
  #                      assays = NULL,
  #                      combine = TRUE)
  #   
  #   # Save the plot
  #   ggplot2::ggsave(filename = paste0("PCA_", j, "_through_", k, ".jpg"),
  #                   plot = last_plot(),
  #                   device = "jpeg",
  #                   path = diagnostics_path,
  #                   scale = 1,
  #                   width = 8.5,
  #                   height = 11,
  #                   units = c("in"),
  #                   dpi = 600,
  #                   limitsize = TRUE,
  #                   bg = "white")
  # }
  
  #**********STEP 10H: DETERMINE THE 'DIMENSIONALITY' FROM ELBOW PLOTS*********#
  
  # # We can quantitatively determine the dimensionality using Elbow plot
  # # We can calculate where the principal components start to elbow by taking the
  # # larger value of:
  # # 1. The point where the principal components only contribute 5% of standard 
  # # deviation and the principal components cumulatively contribute 90% of the 
  # # standard deviation.
  # # 2. The point where the percent change in variation between the consecutive 
  # # PCs is less than 0.1%
  # 
  # # Determine percent of variation associated with each PC
  # pct <- integrated_seurat@reductions$pca@stdev / sum(integrated_seurat@reductions$pca@stdev) * 100
  # # Calculate cumulative percents for each PC
  # cumu <- cumsum(pct)
  # # Determine at which PC the cumulative percent is greater than 90% and 
  # # % variation associated with the PC is less than 5%
  # co1 <- which(cumu > 90 & pct < 5)[1]
  # 
  # # Determine the difference between variation of PC and subsequent PC
  # # We add a 0 to end of diff so that there are equal number of rows in diff, pct & cumu
  # diff <- c(pct[1:length(pct)-1]-pct[2:length(pct)],0) 
  # # Find the last PC which varies from its next PC by more than 0.1%
  # co2 <- max(which(diff > 0.1))
  # 
  # # Minimum of the two calculation
  # pcs <- min(co1, co2+1)
  # cat("\nThe minimum number of PCs to used is ", pcs, "\n")
  # 
  # # Elbow plot to visualize
  # # Create a dataframe with values
  # plot_df <- data.frame(pct = pct,
  #                       cumu = cumu,
  #                       diff = diff,
  #                       rank = 1:length(pct))
  # 
  # ggplot(data = plot_df, aes(x=cumu, y=pct, label=rank, color=rank>pcs)) +
  #   geom_text() +
  #   theme_classic() +       
  #   labs(x = "Cumulative Sum of Std Dev",
  #        y = " Percent Difference in Std Dev",
  #        title = paste0("Elbow plot of PCs")) +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10), 
  #         axis.text.y = element_text(size = 10),                                 
  #         axis.title = element_text(size = 14),                                  
  #         plot.title = element_text(hjust=0.5, size = 16, face="bold")) +       
  #   geom_vline(xintercept = 90, linetype = 2, color = "red") +
  #   geom_hline(yintercept = min(pct[pct > 5]), linetype = 2, color = "red")
  # 
  # ggplot2::ggsave(filename = "Elbow_plot_of_PCA_1_through_50.jpg",
  #                 plot = last_plot(),
  #                 device = "jpeg",
  #                 path = diagnostics_path,
  #                 scale = 1,
  #                 width = 8.5,
  #                 height = 11,
  #                 units = c("in"),
  #                 dpi = 600,
  #                 limitsize = TRUE,
  #                 bg = "white")
  # 
  # # NOTE: We have to use a minimum dimension as determined above. However, using
  # # more PCs is not a bad thing. Clustering will just take more time with more PCs.
  
  return(bad_cells)
}

#******************************************************************************#
#              STEP 14: IDENTIFY CONSERVED MARKERS FOR EACH CLUSTER            #
#******************************************************************************#

# # Set "RNA" assay to be the default assay.
# DefaultAssay(integrated_seurat) <- "RNA"
# 
# # Change active.ident to major clusters we annotated
# Idents(object = integrated_seurat) <- active_res
# 
# # Create function to get conserved markers for any given cluster
# get_conserved <- function(cluster){
#   # min.cells.group = 1 doesnt work in FindConservedMarkers.
#   # So, we manually skip clusters that have less than 3 cells in any sample
#   if (nrow(integrated_seurat@meta.data %>%
#            dplyr::select(Sample, !!rlang::sym(active_res)) %>%
#            dplyr::filter(!!rlang::sym(active_res) == cluster) %>%
#            count(Sample) %>%
#            dplyr::filter(n<3)) == 0){
#     df <- Seurat::FindConservedMarkers(object = integrated_seurat,
#                                        ident.1 = cluster,
#                                        ident.2 = NULL,
#                                        grouping.var = "Sample",
#                                        assay = "RNA",
#                                        slot = "data",
#                                        min.cells.group = 1,
#                                        meta.method = metap::minimump,
#                                        verbose = TRUE,
#                                        features = NULL,
#                                        logfc.threshold = 0.25,
#                                        test.use = "wilcox",
#                                        min.pct = 0.1,
#                                        min.diff.pct = 0.25,
#                                        node = NULL,
#                                        only.pos = TRUE,
#                                        max.cells.per.ident = Inf,
#                                        random.seed = 1,
#                                        latent.vars = NULL,
#                                        min.cells.feature = 3,
#                                        pseudocount.use = 1,
#                                        mean.fxn = NULL,
#                                        fc.name = NULL,
#                                        base = 2,
#                                        return.thresh = 0.01,
#                                        densify = FALSE)
#     
#     if(nrow(df) != 0){
#       df %>%  tibble::rownames_to_column(var="gene") %>%    #add rownames as 1st column named "gene"
#         dplyr::left_join(y = annotations,
#                          by = c("gene" = "SYMBOL")) %>%
#         # merge the dataframe generated by FindConservedMarkers() with unique
#         #genes from "annotations" dataframe columns - gene_name and description
#         cbind(cluster = cluster, .)     # add the cluster variable as 1st column named cluster
#     }
#   }
# }
# 
# # Iterate function across all clusters. We want the output of the map family
# # of functions to be a dataframe with each cluster output bound together by
# # rows, we will use the map_dfr() function
# # map_dfr(inputs_to_function, name_of_function)
# conserved_markers <-
#   map_dfr(rev(levels(integrated_seurat@active.ident)), get_conserved)
# 
# # Save conserved markers
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb = wb1, sheetName = "Conserved_Markers")
# openxlsx::writeData(wb = wb1, sheet = "Conserved_Markers", x = conserved_markers)
# openxlsx::saveWorkbook(wb = wb, file = paste0(seurat_results,"Markers_Conserved.xlsx"), overwrite = TRUE)

#***************(RECOMMENDED): REMOVING BATCH SPECIFIC CLUSTERS**************#

# Sometimes, lots of clusters will be identified which have just 2 cells.
# We can remove them so that UMAPs are not littered with 80+ cluster numbers.

# Sometimes, you will see sample/batch specific clusters where cells of a 
# particular sample dominate. Rather than removing all cells from such a cluster
# we can remove ONLY the cells belonging to dominating sample at all 
# resolutions and re-run entire pipeline for better clustering.

# In example below, we can just remove cells from N6-0, N5-1, N5-11 and N4-12

# Sample  0	            1	            10	          11  	       12
# N1	    40	          38	          240	          0	            130
# N2	    20	          222	          230	          4   	        78
# N3	    632	          350	          2639	        8	            325
# N4	    13	          34	          872	          0	            "5190"
# N5	    399	          "28556"       487	          "6241"	      96
# N6	    "52043"	       562	        1900	        47	          300
# 
# N1	    0.536408743	  0.509588306	  3.218452461	  0	            1.743328416
# N2	    0.109110747	  1.211129296	  1.254773595	  0.021822149	  0.425531915
# N3	    0.854920528	  0.473452824	  3.569834292	  0.010821779	  0.439634765
# N4	    0.098776689	  0.258339032 	6.62563635	  0	            39.43469341
# N5	    0.56627874	  40.52795913	  0.691172296	  8.857507806	  0.136247516
# N6	    71.15045458	  0.768336865	  2.597580149	  0.06425593	  0.410144234
# 
# Total   53147	        29762	        6368	        6300	        6119 

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

#******************************************************************************#
#               OPTIONAL: IDENTIFYING GOOD MARKERS FOR CELL TYPES              # 
#******************************************************************************#

ident_good_markers <- function(){
  # NOTE: You can run this section ONCE to generate a list of good markers and
  # use these markers for ALL SUBSEQUENT data sets. If these markers DO NOT work
  # well for a particular dataset, you can re-run this section using that dataset
  # to generate good markers specific for that data set.
  
  # NOTE: Epithelial Cells and Fibroblasts show vast changes upon BBN condition.
  # So, it might be better to subset the seurat object before finding good markers
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  # Remove unwanted cells and samples before finding good markers
  integrated_seurat1 <- subset(x = integrated_seurat,
                               Condition == "Tumor")
  
  # Create matrix similar to read data for each cluster to identify closely related clusters
  cor_meta_data <- integrated_seurat1@meta.data %>%
    dplyr::distinct(integrated_snn_res.1.4, .keep_all = TRUE) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("integrated_snn_res.1.4")
  
  clusters <- as.vector(unlist(integrated_seurat1@meta.data %>% dplyr::distinct(integrated_snn_res.1.4)))
  cor_read_data <- matrix(NA, nrow = nrow(integrated_seurat1@assays$RNA@counts), ncol = length(clusters))
  rownames(cor_read_data) <- rownames(integrated_seurat1@assays$RNA@counts)
  colnames(cor_read_data) <- clusters
  for(i in 1:length(clusters)){
    
    # Create a list of cells for each sample
    cor_cells_subset <- rownames(integrated_seurat1@meta.data %>% dplyr::filter(integrated_snn_res.1.4 == as.numeric(clusters[i])))
    
    # Get the raw counts of these cells. counts are stored in sparse matrix format
    # So, use data.frame to convert "." in sparse matrix to "0"
    cor_subset <- data.frame(integrated_seurat1@assays$RNA@counts[,cor_cells_subset])
    cor_read_data[,i]  <- rowSums(cor_subset)
  }
  
  # Remove unwanted genes
  unwanted_genes <- rownames(cor_read_data)[grep(pattern = "^RP[SL]|RIK$|^MT-|^GM[0-9.]+$", x= rownames(cor_read_data), ignore.case=TRUE)]
  keep_genes <- setdiff(rownames(cor_read_data), unwanted_genes)  
  cor_read_data <- cor_read_data[keep_genes,]
  
  # Create DESeq2 object so we can use DESeq2::rlog() or DESeq2::vst() on it
  cor_dds <- DESeqDataSetFromMatrix(countData = cor_read_data,
                                    colData = cor_meta_data,
                                    design = ~1)
  # Use vst or rld
  vst <- DESeq2::vst(cor_dds)
  head(assay(vst))
  sampleDists <- dist(t(assay(vst)))
  sampleDistMatrix <- as.matrix(sampleDists)
  pheatmap::pheatmap(mat = sampleDistMatrix,
                     scale = "none",
                     cutree_rows = 20,
                     cutree_cols = 20,
                     cluster_rows = TRUE,   #cluster the rows
                     cluster_cols = TRUE,   #cluster the columns
                     fontsize = 8, 
                     fontsize_row = 8, 
                     fontsize_col = 8,
                     angle_col = c("270", "0", "45", "90", "315"),
                     fontsize_number = 0.8*fontsize,
                     filename = paste0(seurat_results, "Cluster_Correlation", "_", celltype, ".jpg"))
  
  # Use the UMAP at res 1.4 as well as above correlation heatmap to group clusters
  res <- "integrated_snn_res.1.4"
  
  # BBN_C57B6 scRNASeq res 1.4 (only BBN samples)
  # Refer ppt to understand how the clusters were grouped
  clusters <- list("cluster1" = c(13,16,17,23,29,35),     #Fibro_1
                   "cluster2" = c(5,6,26,30,33,42),       #Fibro_2
                   "cluster3" = c(36),                    #Myocytes
                   "cluster4" = c(0,8,11,14,15,18,25,34,39), #Epithelial
                   "cluster5" = c(4,7,12,21),             #MDSCs
                   "cluster6" = c(9,10,20),               #Macrophages
                   "cluster7" = c(32,37),                 #DCs
                   "cluster8" = c(19,22,24,27,28),        #T cells
                   "cluster9" = c(31),                    #NK cells
                   "cluster10" = c(2,48),                 #B cells
                   "cluster11" = c(1,46),                 #Endo
                   "cluster12" = c(38),                   #Lymp. Endo
                   "cluster13" = c(43),                   #Neurons??
                   "mixed" = c(3,40,41,50,51,52,53,54,44,45,47,49)) #contaminants
  
  # BBN_C57B6 scRNASeq Epithelial res 1.4 (only BBN samples)
  # Refer ppt to understand how the clusters were grouped
  clusters <- list("cluster1" = c(4,18),     
                   "cluster2" = c(1,8,9,20,21),      
                   "cluster3" = c(3,12,19,22), 
                   "cluster4" = c(7,23),             
                   "cluster5" = c(2,10),               
                   "cluster6" = c(0,15),                 
                   "cluster7" = c(5,11,17),        
                   "cluster8" = c(6,16),                    
                   "cluster9" = c(13),                 
                   "cluster10" = c(14),                 
                   "mixed" = c(24,25,26,27,28,30))
  
  # BBN_C57B6 scRNASeq Fibroblast res 1.4 (only BBN samples)
  # Refer ppt to understand how the clusters were grouped
  clusters <- list("cluster1" = c(0,4,5),     
                   "cluster2" = c(3,6,12,14,18),      
                   "cluster3" = c(7,8,9,13,15), 
                   "cluster4" = c(2,10,17),             
                   "cluster5" = c(11),               
                   "cluster6" = c(1,16),                 
                   "cluster7" = c(20),
                   "mixed" = c(19))
  
  # BBN_C57B6 scRNASeq Myeloid res 1.4 (only BBN samples)
  # Refer ppt to understand how the clusters were grouped
  clusters <- list("cluster1" = c(2,4),     
                   "cluster2" = c(0,1,8,11),      
                   "cluster3" = c(5,6,9,10,19,20), 
                   "cluster4" = c(3,7,12,21,24),             
                   "cluster5" = c(13,14,15,16,22),               
                   "cluster6" = c(17), 
                   "cluster7" = c(18), 
                   "cluster8" = c(25),
                   "cluster9" = c(26),
                   "mixed" = c(23,27,28,29))
  
  
  # Make sure you have assigned all clusters to one of the cell types
  # We subtract 1 because (...find explanation..)
  list_1 <- sort(unique(as.numeric(unlist(integrated_seurat1@meta.data %>% 
                                            dplyr::select(all_of(res)), use.names=FALSE))))-1
  list_2 <- sort(unlist(clusters, use.names=FALSE))
  
  # Proceed with annotation ONLY if all clusters have been renamed
  if (identical(list_1, list_2)){
    "All Clusters have been annotated"
    
    # Extract metadata from Seurat object, assign appropriate resolution to
    # seurat_clusters column and add cell_type column
    data <- integrated_seurat1@meta.data %>% 
      dplyr::mutate(seurat_clusters = !!rlang::sym(res),
                    cell_type = NA)
    
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
    integrated_seurat1@meta.data <- data
  } else {
    "Some clusters have not been annotated yet"
  }
  
  # STEP 2: Run FindAllMarkers() on these major clusters
  
  # Set "RNA" assay to be the default assay.
  DefaultAssay(integrated_seurat1) <- "RNA"
  
  # Change active.ident to major clusters we annotated
  Idents(object = integrated_seurat1) <- "cell_type"
  
  #Remove unwanted cells and samples before finding good markers
  integrated_seurat1 <- subset(x = integrated_seurat1,
                               subset = (cell_type == "mixed"), 
                               invert=TRUE)
  
  # Get all markers for any given cluster
  all_markers <- Seurat::FindAllMarkers(object=integrated_seurat1,
                                        assay = "RNA",
                                        features = NULL,
                                        logfc.threshold = 0.25,
                                        test.use = "wilcox",
                                        slot = "data",
                                        min.pct = 0.1,
                                        min.diff.pct = 0.1,
                                        node = NULL,
                                        verbose = TRUE,
                                        only.pos = TRUE,
                                        max.cells.per.ident = Inf,
                                        random.seed = 1,
                                        latent.vars = NULL,
                                        min.cells.feature = 3,
                                        min.cells.group = 1,
                                        pseudocount.use = 1,
                                        mean.fxn = NULL,
                                        fc.name = NULL,
                                        base = 2,
                                        return.thresh = 0.01,
                                        densify = FALSE)
  
  all_markers <- all_markers %>% 
    dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio = pct.1/pct.2) %>%
    dplyr::left_join(y = unique(annotations[, c("SYMBOL", "chr", "description")]), by = c("gene" = "SYMBOL")) %>%
    dplyr::relocate(cluster, gene, chr, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, description)
  
  # STEP 3: Find top markers for each major cluster
  # (i) statistically significant (p_val_adj < 0.05)
  # (ii) low expression in other cell types (high ratio)
  # (iii) expressed at decent levels (avg_log2FC >= 0.58)
  # NOTE: We are not using pct.1 as a good marker may not be expressed in most of
  # cells in the major cluster as we have merged many smaller clusters together
  top_markers <- all_markers %>%
    dplyr::filter(avg_log2FC >= 0.58 & p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(ratio), desc(avg_log2FC)) %>%
    dplyr::slice_head(n=20) %>%
    ungroup()
  
  # Save all the markers
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "All_Markers")
  openxlsx::writeData(wb = wb, sheet = "All_Markers", x = all_markers)
  openxlsx::addWorksheet(wb = wb, sheetName = "Top_Markers")
  openxlsx::writeData(wb = wb, sheet = "Top_Markers", x = top_markers)
  openxlsx::saveWorkbook(wb = wb, file = paste0(seurat_results, "Markers_for_major_clusters_", celltype, ".xlsx"), overwrite = TRUE)
  
  # Step 4: Manually rename clusters with appropriate cell types using info 
  # below in the xlsx file and use it for future annotations. Save this info
  # to "KMS_Markers.xlsx" in sheet named "All Markers"
  
  # Fibroblasts express several types of collagens (Col*)
  # Epithelial cells express several types of keratins (Krt*)
  # T cells express CD3 genes (Cd3d, Cd3e, Cd3g)
  # B cells express CD79 genes (Cd79a, Cd79b)
  # Macrophages express C1q genes (C1qa, C1qb, C1qc)
  # Mast Cells express S100A8, S100A9, Retnlg at high levels
  # Myocytes express myosin genes (Myh11, Myl9)
  # Endothelial cells express endothelial proteins (Esam, Ecsr, Emcn, Nos3)
  # Neurons express several potassium gated channels (Kcna1, Kcna2, Kcna6)
  
  # Look at UMAPs plotted from FeaturePlot() below and remove bad markers.
  
  features <- read.xlsx(paste0(seurat_results, "Markers_for_major_clusters_", celltype, ".xlsx"), sheet="Top_Markers")
  features <- feature$gene
  features <- intersect(features, rownames(integrated_seurat@assays$RNA@data))
  length(features)
  
  for(i in 1:ceiling(length(features)/10)){
    j <- 10*i-9
    k <- 10*i
    Seurat::FeaturePlot(object = plot_object,  #integrated_seurat,
                        features = features[j:k],
                        pt.size = 0.4,
                        order = TRUE,
                        min.cutoff = 'q10',
                        reduction = "umap",
                        cols =  c("grey", viridis(n = 10, option = "C", direction = -1)),
                        ncol =5,
                        label = FALSE,
                        combine = TRUE,
                        raster = FALSE)
    
    ggplot2::ggsave(filename = paste0(i,".jpg"),
                    plot = last_plot(),
                    device = "jpeg",
                    #path = diagnostics_path,
                    scale = 1,
                    width = 22,
                    height = 11,
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = "white")
  }
}

#******************************************************************************#
#                                  RUN DESEQ2                                  #
#******************************************************************************#

# This function runs DESeq2 and plots volcano plot, heatmap, PCA plot and
# correlation plot.
# meta_data MUST have "Sample" column with sample names without any duplication
# meta_data MUST have "Batch" column with batch info
# read_data MUST have "SYMBOL" column with gene names without any duplication
# Several global variables need to be defined. See RNASeq.R

analyze_DESeq2 <- function(meta_data, read_data, file_suffix){
  
  #****************************************************************************#
  #                           STEP 1: IMPORT META DATA                         #
  #****************************************************************************#
  
  # SUBSET METADATA (VERY VERY IMPORTANT, OFTEN OVERLOOKED ASPECT OF ANALYSIS)
  # DESeq2 calculates sizefactors using all the samples in meta_data. So, final
  # results are heavily dependent on sizefactors. For instance, if you have 
  # 5 replicates of Vehicle treated samples, 5 replicates of Enzalutamide 
  # treated samples and 5 replicates of R1881 treated samples, and you want to 
  # find DEGs between R1881 vs vehicle, remove the Enzalutamide treated samples 
  # from meta_data & read_data. Since, we have to apply batch correction on the
  # whole dataset, we keep all samples for now but we remove them before 
  # differential expression analysis.
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#if-i-have-multiple-groups-should-i-run-all-together-or-split-into-pairs-of-groups
  
  meta_data <- meta_data %>%
    dplyr::mutate(id = dplyr::if_else(get(Variable) %in% c(Comparisons$Target, Comparisons$Reference) & 
                                        get(Variable2) %in% Variable2_value, 
                                      paste0(get(Variable), "_", get(Variable2)), 
                                      "ignore")) %>%
    dplyr::mutate(id = dplyr::if_else(is.na(id), "ignore", id)) %>%
    dplyr::filter(id != "ignore") %>% 
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var="Sample")
  
  #****************************************************************************#
  #                          STEP 2: IMPORT READ DATA                          #                                           
  #****************************************************************************#
  
  read_data <- read_data %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var="SYMBOL") %>%
    dplyr::select(intersect(colnames(.), rownames(meta_data))) %>%
    replace(is.na(.), 0)
  
  #****************************************************************************#
  #              STEP 3: PREPARE META DATA and READ DATA FOR DESEQ2            #
  #****************************************************************************#
  
  # DESeq2 automatically removes genes that have 0 counts. However, we have to
  # make sure there are no samples which have 0 counts for all genes. For eg,
  # sample MB1 has 0 counts for all genes in cluster 35. So, we need to manually
  # remove this sample before DESeq2 analysis. Else, the geometric mean will be 
  # zero for all genes and DESeq2 will halt execution.
  
  if (TRUE %in% (colSums(read_data)==0)){
    read_data <- read_data[, -c(which(colSums(read_data)==0))]
    meta_data <- meta_data[-c(which(colSums(read_data)==0)),]
  }
  
  # Rearrange read_data to match meta_data
  read_data <- read_data[,rownames(meta_data)]
  
  # Remove samples that have "NA" for variable being analyzed in metadata.
  NA_sample_list <- which(is.na(meta_data[Variable]))
  if (length(NA_sample_list)!=0){
    meta_data <- meta_data[-c(NA_sample_list),]
    read_data <- read_data[,-c(NA_sample_list)]
  }
  
  # Remove any column named "sizeFactor" in meta_data as it interferes with DESeq2()
  if ("sizeFactor" %in% colnames(meta_data)){
    meta_data <- meta_data %>% dplyr::select(!sizeFactor)
  }
  
  # Check if all columns of meta table are factors. Else, convert to factors.
  print("meta_data before factor conversion")
  str(meta_data)
  for (i in 1:ncol(meta_data)){
    meta_data[,i] <- as.factor(meta_data[,i])
  }
  print("meta_data after factor conversion")
  str(meta_data)
  
  # Check (i) read_data & meta_data are dataframes
  # (ii) column names of read_data (i.e.sample name) are PRESENT in row names of meta_data (i.e. sample name)
  # (iii) column names of read_data (i.e.sample name) are PRESENT in same order as row names of meta_data (i.e. sample name)
  if (class(read_data) == "data.frame" & class(meta_data) == "data.frame" & 
      base::all(colnames(read_data) %in% rownames(meta_data)) &
      base::all(colnames(read_data) == rownames(meta_data))){
    print("ALL OK")
  } else if (class(read_data) != "data.frame"){
    print("read_data is NOT a dataframe")
  } else if (class(meta_data) != "data.frame"){
    print("meta_data is NOT a dataframe")
  } else if (!(all(colnames(read_data) %in% rownames(meta_data)))){
    print("Missing samples")
  } else if (all(colnames(read_data) != rownames(meta_data))){
    print("Varying sample order")
  } else{}
  
  #****************************************************************************#
  #                      STEP 4: PERFORM BATCH CORRECTION                      #
  #****************************************************************************#
  
  if (length(unique(as.numeric(meta_data$Batch))) > 1){
    # Adjust for batch effects using Combat-Seq
    corrected_read_data <- sva::ComBat_seq(counts = as.matrix(read_data), 
                                           batch = as.numeric(meta_data$Batch), 
                                           group = as.numeric(as.factor(meta_data$id)))
  } else{
    corrected_read_data <- read_data
  }
  
  #****************************************************************************#
  #                     STEP 5: CALCULATE NORMALIZED COUNTS                    #
  #****************************************************************************#
  
  # NOTE: sizeFactors MUST be estimated first before normalization.
  # sizeFactors will differ based on samples included i.e. all samples vs subset
  # of samples. This will in turn affect normalized counts.
  # sizeFactors doesnt depend on the design formula. So, design ~1 or design ~id
  # will give same normalized values.
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  
  normalized_counts <- counts(dds, normalized=TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ENSEMBL_ID")
  
  if (length(intersect(annotations$ENSEMBL_ID, normalized_counts$ENSEMBL_ID)) > length(intersect(annotations$SYMBOL, normalized_counts$ENSEMBL_ID))){  
    normalized_counts <- annotations %>%
      dplyr::right_join(normalized_counts, by = c("ENSEMBL_ID"="ENSEMBL_ID"), multiple = "all") %>%
      dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                              is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                              TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("ENSEMBL_ID", .keep_all = TRUE)
  } else{
    normalized_counts1 <- annotations %>%
      dplyr::right_join(normalized_counts, by = c("SYMBOL"="ENSEMBL_ID"), multiple = "all") %>%
      dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                              is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                              TRUE ~ SYMBOL)) %>% 
      dplyr::distinct_at("SYMBOL", .keep_all = TRUE)
  }
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "normalized")
  openxlsx::writeData(wb, sheet = "normalized", x = normalized_counts, rowNames = FALSE)
  
  openxlsx::addWorksheet(wb, sheetName = "raw")
  openxlsx::writeData(wb, sheet = "raw", x = corrected_read_data, rowNames = FALSE)
  openxlsx::saveWorkbook(wb,
                         file = paste0(results_path, "Results_", file_suffix, "_", Variable2_value, "_Normalized_Counts.xlsx"),
                         overwrite = TRUE)
  
  #**************************************************************************#
  #                          STEP 6: PERFORM DESEQ2                          #
  #**************************************************************************#
  
  meta <- meta_data
  read <- corrected_read_data
  
  # Perform DEG analysis for each comparison with the dataset
  for (n in 1:length(Comparisons$Target)){
    
    # Subset meta_data
    meta_data <- meta %>% dplyr::filter(get(Variable) %in% c(Comparisons$Target[n], Comparisons$Reference[n]))
    
    # Subset read_data
    corrected_read_data <- read %>% dplyr::select(rownames(meta_data))
    
    #**************************************************************************#
    #                          STEP 5: PERFORM DESEQ2                          #
    #**************************************************************************#
    
    # Create DESeq2Dataset object, estimate Size factors and view it
    # NOTE: sizefactors MUST be calculated ONLY using samples being compared for 
    # differential expression. So, make sure read_data and meta_data ONLY have
    # samples being compared for differential expression.
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ id)
    
    # Study effect of Treatment ignoring effect of Sex
    #dds <- DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ Treatment)
    
    # Study effect of Treatment after removing effect of Sex (assumes equal effect of Sex on different treatments)
    #dds <- DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ Sex + Treatment)
    
    # Study effect of Treatment after removing effect of Sex (assumes different effect of Sex on different treatments)
    #dds <- DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ Sex + Treatment + Sex:Treatment)
    
    # Determine whether to use "local" or "parametric" fit
    dds_para <- DESeq2::DESeq(object = dds, 
                              test = "Wald",
                              fitType = "parametric",
                              betaPrior = FALSE,
                              minReplicatesForReplace = 7)
    dds_local <- DESeq2::DESeq(object = dds, 
                               test = "Wald",
                               fitType = "local",
                               betaPrior = FALSE,
                               minReplicatesForReplace = 7)
    residual_para <- mcols(dds_para)$dispGeneEst - mcols(dds_para)$dispFit
    residual_local <- mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit
    if (median(residual_para^2, na.rm=TRUE) <= median(residual_local^2, na.rm=TRUE)){
      fit <- "parametric"
    } else{
      fit <- "local"
    }
    
    # Run analysis
    # betaPrior: default=FALSE, shrunken LFCs are obtained later using lfcShrink
    # fitType: "local" may be better sometimes
    dds <- DESeq2::DESeq(object = dds, 
                         test = "Wald",
                         fitType = fit,
                         betaPrior = FALSE,
                         minReplicatesForReplace = 7)
    
    # Define contrast and coeff
    DE_levels <- as.vector(unique(meta_data$id))
    contrast <- c("id", 
                  DE_levels[grepl(pattern = Comparisons$Target[n], x=DE_levels)],
                  DE_levels[grepl(pattern = Comparisons$Reference[n], x=DE_levels)])
    
    coeff <- paste(contrast[1], contrast[2], "vs", contrast[3], sep='_')
    cat("Coeff is ", coeff, "\n")
    
    # Calculate results
    res <- DESeq2::results(object = dds,
                           contrast = contrast,
                           lfcThreshold = lfc.cutoff,
                           altHypothesis = "greaterAbs",
                           cooksCutoff = TRUE,
                           independentFiltering = TRUE,
                           alpha = padj.cutoff,
                           pAdjustMethod = "BH")
    
    # Perform lfcshrinkage to account for variability between replicates
    # For ashr, if res is provided, then coef and contrast are ignored.
    # lfcshrinkage will not change the number of DEGs and affects only logFC
    res <- DESeq2::lfcShrink(dds = dds, res = res, type = "ashr")
    
    # Summarize results
    summary(res)
    
    # Export results as dataframe and replace ensembl IDs with Gene names
    results <- res %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ENSEMBL_ID") 
    
    if (length(intersect(annotations$ENSEMBL_ID, results$ENSEMBL_ID)) > length(intersect(annotations$SYMBOL, results$ENSEMBL_ID))){  
      results <- annotations %>%
        dplyr::right_join(results, by = c("ENSEMBL_ID"="ENSEMBL_ID"), multiple = "all") %>%
        dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
        dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                                is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                                TRUE ~ SYMBOL)) %>%
        dplyr::distinct_at("ENSEMBL_ID", .keep_all = TRUE)
    } else{
      results <- annotations %>%
        dplyr::right_join(results, by = c("SYMBOL"="ENSEMBL_ID"), multiple = "all") %>%
        dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
        dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                                is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                                TRUE ~ SYMBOL)) %>%
        dplyr::distinct_at("SYMBOL", .keep_all = TRUE)
    }
    
    # Save the results
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "DEGs")
    openxlsx::writeData(wb, sheet = "DEGs", x = results, rowNames = FALSE)
    openxlsx::saveWorkbook(wb, 
                           file = paste0(results_path, "Results_", file_suffix, "_", Variable2_value, "_", coeff, "_DEGs.xlsx"), 
                           overwrite = TRUE)
    
    #**************************************************************************#
    #             STEP 6: CREATE HEATMAPS, VOLCANO PLOTS, PCA PLOTS            #
    #**************************************************************************#
    
    if (heatmap_plot == TRUE){
      
      # Refer https://www.biostars.org/p/339065/ where the authors of DESeq2 recommend
      # using vst or rld rather than log1p transformation while plotting DESeq2 
      # normalized counts in heatmap.
      # NOTE: We scale the vst values and replace all NAs with 0 before feeding it to 
      # pheatmap(). If we perform the scaling in pheatmap, the NAs generated during
      # scaling will give error.
      # NOTE: While use of log1p transformation is NOT recommended, I notice only 
      # minor changes in color
      
      # Define any filename you want added to final file
      file_suffix <- paste0("_", Variable2_value, "_", coeff)
      
      # (I) Read metadata
      metadata <- meta_data %>%
        tibble::rownames_to_column("Sample")
      
      # (II) Define genes to plot
      plot_genes <- results %>%
        dplyr::filter(padj < 0.05) %>% 
        dplyr::select(SYMBOL) %>%
        unlist(use.names = FALSE)
      
      # (III) Genes to display in heatmap
      disp_genes <- c()
      
      # (IV) Read expr data
      #normalized_counts <- normalized_counts
      normalized_counts <- assay(DESeq2::vst(dds, blind=FALSE)) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("ENSEMBL_ID") %>%
        dplyr::left_join(annotations, by = c("ENSEMBL_ID"="ENSEMBL_ID"), multiple = "all") %>%
        dplyr::filter(nchar(SYMBOL) > 0) %>%
        dplyr::select(SYMBOL, everything(), -c(ID, CHR, DESCRIPTION, START, END, STRAND))
      colnames(normalized_counts)[1] <- "SYMBOL"
      
      # Run the function
      plot_heatmap(normalized_counts, metadata, plot_genes, disp_genes, file_suffix)
    }
    
    if (volcano_plot ==TRUE){
      
      # Define any filename you want added to final file
      file_suffix <- paste0("_", Variable2_value, "_", coeff)
      
      # (I) Import expression data with log2FC and pval
      volcano_df <- results %>%
        dplyr::rename(log2FC = log2FoldChange, Gene = SYMBOL)
      
      # (II) Define or import metadata
      # NOTE: Metadata is a dataframe with "Sample" and "Condition" columns
      metadata <- meta_data %>%
        tibble::rownames_to_column("Sample")
      
      # (III) Define any genes you want to mark in volcano plot
      disp_genes <- c()
      
      # Make volcano plots
      plot_volcano(volcano_df, disp_genes,file_suffix)
    }
    
    if (PCA_plot == TRUE){
      
      # Define any filename you want added to final file
      file_suffix <- paste0("_", Variable2_value, "_", coeff)
      
      # Plot PCA
      rld_mat <- assay(DESeq2::rlog(dds, blind=FALSE))
      pca <- prcomp(t(rld_mat))
      df <- cbind(meta_data, pca$x)      # Create data frame with metadata, PC1 & PC2 values for ggplot
      ggplot(df) + geom_point(aes(x=PC1, y=PC2, color = Condition))
      ggplot2::ggsave(filename =  paste0("PCA_Plot", file_suffix, ".pdf"),
                      plot = last_plot(),
                      device = "pdf",
                      path = results_path,
                      scale = 1,
                      #width = 8.5,
                      #height = 11,
                      units = c("in"),	 
                      dpi = 300,
                      limitsize = TRUE,
                      bg = "white")
    }
    
    if (cor_plot==TRUE){
      
      #Determine number of cuts in heatmap
      cuts <- 20
      
      # Use vst as well as rld
      vst <- DESeq2::vst(dds)
      vst <- vst[rownames(vst) %in% rownames(results %>% dplyr::filter(abs(log2FoldChange) >= 0.58 & padj < 0.05 & padj > 0 & !is.na(padj)))]
      head(assay(vst))
      sampleDists <- dist(assay(vst))
      sampleDistMatrix <- as.matrix(sampleDists)
      rownames(sampleDistMatrix) <- dplyr::left_join(data.frame("ENSEMBL_ID" = rownames(sampleDistMatrix)), annotations, by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>% 
        dplyr::select(SYMBOL) %>% 
        unlist(use.names = FALSE)
      colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
      
      # cluster and re-order rows
      rowclust <- hclust(dist(sampleDistMatrix))
      
      # cluster and re-order columns
      colclust <- hclust(dist(t(sampleDistMatrix)))
      reordered <- sampleDistMatrix[rowclust$order,colclust$order]
      
      # Save batch corrected normalized counts for entire dataset
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, sheetName = "Correlation")
      openxlsx::writeData(wb, sheet = "Correlation", x = reordered , rowNames = TRUE)
      openxlsx::saveWorkbook(wb,
                             file = paste0(results_path, "Results_", Variable2_value, "_Correlation_Matrix.xlsx"),
                             overwrite = TRUE)
      
      pheatmap::pheatmap(mat = reordered,
                         scale = "none",
                         cellwidth = 3,
                         cellheight = 2,
                         cutree_rows = cuts,
                         cutree_cols = cuts,
                         cluster_rows = TRUE,   #cluster the rows
                         cluster_cols = TRUE,   #cluster the columns
                         fontsize = 8, 
                         fontsize_row = 8, 
                         fontsize_col = 8,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         angle_col = c("270", "0", "45", "90", "315"),
                         fontsize_number = 0.8*fontsize,
                         width = 40,
                         height = 40,
                         filename = paste0(results_path, "Correlation_vst.jpg"))
    }
  }
}

#******************************************************************************#
#         PREPARE READ AND META DATA FROM SEURAT OBJECT AND RUN DESEQ2         #
#******************************************************************************#

# This function reads a seurat object and generates read_data, meta_data, file_suffix
# for use by analyze_DESeq2(). It automatically calls analyze_DESeq2() too.

prep_DESeq2 <- function(celltype){ 
  
  cat("\n**********************", celltype, "Analysis**********************\n")
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  #******************************IMPORT META DATA******************************#
  
  integrated_seurat <- subset(x = integrated_seurat,
                              subset = cell_class %in% c("Mixed", "Unclassified"),
                              invert = TRUE)
  
  # Perform analysis on cell types like T cell, B cell rather than Myeloid, Lymphoid
  
  subtypes <- integrated_seurat@meta.data %>% 
    dplyr::select(cell_type) %>% 
    unlist(use.names=FALSE) %>%
    unique()
  
  for (subtype in subtypes){
    
    subset_seurat <- subset(x = integrated_seurat,
                            subset = cell_type == subtype)
    
    # Subset metadata
    meta_data <- subset_seurat@meta.data %>%
      dplyr::distinct(Sample, .keep_all = TRUE) %>%
      dplyr::mutate(Batch = 1)
    
    #******************************IMPORT READ DATA******************************#                                           
    
    # The read data will have "the reads of all cells belonging to a single 
    # sample" merged together in each column. First, create a list of samples
    samples <- subset_seurat@meta.data %>% 
      dplyr::select(Sample) %>% 
      unlist(., use.names=FALSE) %>% 
      unique()
    
    # Second, create an empty dataframe with rows=genes and columns=samples
    read_data <- data.frame(matrix(NA, nrow = nrow(subset_seurat@assays$RNA@counts), ncol = nrow(meta_data)))
    rownames(read_data) <- rownames(subset_seurat@assays$RNA@counts)
    colnames(read_data) <- samples
    
    # Thirdly, we will add row-wise, the counts of each gene for each sample
    for(i in samples){
      
      # Create a list of cells for each sample
      cells_subset <- rownames(subset_seurat@meta.data %>% dplyr::filter(Sample == i))
      
      # Use data.frame to convert "." in sparse matrix to "0"
      subset <- data.frame(subset_seurat@assays$RNA@counts[,cells_subset])
      read_data[,i]  <- rowSums(subset)
    }
    
    read_data <- read_data %>% 
      tibble::rownames_to_column("SYMBOL")
    
    file_suffix <- subtype
    analyze_DESeq2(meta_data, read_data, file_suffix)
  }
}

#*********************************CLEAR MEMORY*********************************#

get_homolog <- function(){
  mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  convert_mouse_to_human <- function(gene_list){
    
    output = c()
    
    for(gene in gene_list){
      class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
      if(!identical(class_key, integer(0)) ){
        human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
        for(human_gene in human_genes){
          output = append(output,human_gene)
        }
      }
    }
    
    return (output)
  }
}
#******************************************************************************#