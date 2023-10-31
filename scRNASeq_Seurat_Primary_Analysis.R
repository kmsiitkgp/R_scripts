#!/usr/bin/env Rscript

# NOTE: All variables and functions are defined within the file below

source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

#******************************************************************************#
#                       STEP 1: SETUP THE SEURAT OBJECT                        #
#******************************************************************************#

#************************IMPORTING DATA FROM h5AD FILE*************************#

if (h5ad_file == "yes"){
  
  # Load h5ad (useful if analyzing collaborator data in h5ad format)
  SeuratDisk::Convert(source = paste0(seurat_results, proj, ".h5ad"),
                      dest = "h5seurat",
                      assay="RNA",
                      overwrite = FALSE)
  
  raw_seurat <- SeuratDisk::LoadH5Seurat(file = paste0(seurat_results, proj, ".h5seurat"))
  
  # View the contents of the seurat object
  dplyr::glimpse(raw_seurat)
}

#***********IMPORTING DATA FROM OUTPUT FOLDERS OF CELLRANGER/CITESEQ***********#

if (demultiplexed == "yes"){
  
  # Create a list of samples which will be added to each barcode.
  # Since folder names correspond to sample name, we just use list.files()
  samples <- list.files(path = demux_results)
  
  for (i in samples){
    sample.seurat <- readRDS(file = paste0(demux_results, i))
    assign(i, sample.seurat)
  }
} 

if (demultiplexed == "no"){
  
  # Create a list of samples which will be added to each barcode.
  # Since folder names correspond to sample name, we just use list.files()
  samples <- list.files(path = feature_matrix_path)
  
  # Loop through each of the individual folders in parent directory & import data
  for(i in samples){
    
    # Read the data files from each sample folder
    # NOTE: gene.column=1 imports Ensembl ids from features.tsv. DO NOT DO THIS. 
    # Read gene symbols from features.tsv as we will calculate mitoratio, 
    # riboratio etc using gene names
    sample.dgCMatrix <- Seurat::Read10X(data.dir = paste0(feature_matrix_path, i),
                                        gene.column = 2,  
                                        cell.column = 1,
                                        unique.features = TRUE,
                                        strip.suffix = FALSE)
    
    # Create a seurat object for each dgCMatrix object
    # NOTE: Set min.features=100 to remove low quality cells & reduce object size
    sample.seurat <- SeuratObject::CreateSeuratObject(counts = sample.dgCMatrix,
                                                      project = i,
                                                      assay = "RNA",
                                                      names.field = 1,
                                                      names.delim = "_",
                                                      meta.data = NULL,
                                                      min.cells = 0,
                                                      min.features = 100,
                                                      row.names = NULL)
    
    # Assign the seurat object to its corresponding variable
    assign(i, sample.seurat)
    
    # Explore the meta.data slot
    cat("\nFirst few rows of ", i, "\n")
    print(head(sample.seurat@meta.data))
    cat("\nLast few rows of ", i, "\n")
    print(tail(sample.seurat@meta.data))
  }
}

# Create a merged Seurat object.
# NOTE: We are going to do the same QC on each sample. So, we merge the 
# individual seurat objects into a single seurat object. However, samples will
# have same barcodes. To keep track of cell identities (i.e.barcodes) coming 
# from each sample after merging, we add a prefix (i.e. sample name) to each 
# barcode using add.cell.ids
raw_seurat <- base::merge(x = get(samples[1]),
                          y = lapply(samples[2:length(samples)], get),
                          add.cell.ids = samples,
                          merge.data = FALSE)

#******************************************************************************#
#                           STEP 2: QUALITY CONTROL                            #
#******************************************************************************#

#**********************STEP 2A: CALCULATE ALL QC METRICS***********************#

# Compute percent mito percent
raw_seurat <- Seurat::PercentageFeatureSet(object = raw_seurat,
                                           pattern = dplyr::if_else(species=="Homo sapiens", "MT-", "mt-"),
                                           features = NULL,
                                           col.name = "MitoPercent",
                                           assay = NULL)

# Compute percent ribo percent
raw_seurat <- Seurat::PercentageFeatureSet(object = raw_seurat,
                                           pattern = dplyr::if_else(species=="Homo sapiens", "^RP[SL]", "^Rp[sl]"),
                                           features = NULL,
                                           col.name = "RiboPercent",
                                           assay = NULL)

# Extract metadata
raw_metadata <- raw_seurat@meta.data

# Rename columns to be more intuitive and add the additional QC metrics:
# (i)     Cell      : Unique identifiers corresponding to each cell = barcodes = row names of metadata dataframe
# (ii)    Sample    : sample name
# (iii)   nUMIs     : number of transcripts per cell
# (iv)    nGenes    : number of genes per cell
# (v)     nHTO_UMIs : number of HTO reads per cell
# (vi)    nHTOs     : number of HTOs types per cell
# (Vii)   MitoRatio : MitoPercent/100
# (viii)	RiboRatio : RiboPercent/100  
# (ix)    Novelty   : log ratio of genes per UMI
# NOTE: B2 was demuxed with HTODemux. So, we use HTO_classification ONLY for B2.
raw_metadata <- raw_metadata %>% 
  dplyr::transmute(Cell = rownames(raw_metadata),
                   Sample = orig.ident,
                   #Unique_ID = dplyr::if_else(Sample == "B2", paste0(Sample, "_", HTO_classification), paste0(Sample, "_",  MULTI_classification)),
                   nUMIs = nCount_RNA,
                   nGenes = nFeature_RNA,
                   #nHTO_UMIs = nCount_HTO,
                   #nHTOs = nFeature_HTO,
                   #HTO_tag = MULTI_classification,
                   MitoRatio = MitoPercent/100,
                   RiboRatio = RiboPercent/100,
                   Novelty = log10(nGenes)/log10(nUMIs))

# Import other meta data associated with data set. 
# NOTE: This xslx file should "ONLY" have necessary metadata that you want to 
# add like Sex, Treatment, etc along with 
extra_metadata <-  openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, metafile),
                                       colNames = TRUE)

# Merge imported metadata with existing metadata
# NOTE: left_join will remove rownames. So, add rownames before replacing 
# metadata in Seurat object
if (demultiplexed == "yes"){
  # Sample was initially assigned from orig.ident which corresponds to batch.
  # So, we import sample information based on hashtag identified from metadata
  raw_metadata <- raw_metadata %>%
    dplyr::select(everything(), -c(Sample)) %>%
    dplyr::left_join(extra_metadata, by=("HTO_tag"="HTO_tag")) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::column_to_rownames(var = "index")
} else{
  raw_metadata <- raw_metadata %>% 
    dplyr::left_join(extra_metadata, by=("Sample"="Sample")) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::column_to_rownames(var = "index")
}

# Replace the metadata in raw Seurat object
raw_seurat@meta.data <- raw_metadata

#******************************STEP 2B: PERFORM QC*****************************#

# # Calculate number of cells that remain at different QC settings 
# summary_df <- data.frame(Sample = "NA")
# cutoffs <- list(gene_cutoff  = c(250, 500, 1000, 5000),
#                 umi_cutoff = c(500, 1000, 10000, 100000),
#                 mito_cutoff = c(0.05, 0.1, 0.2, 0.3),
#                 novelty_cutoff = c(0.8, 0.85, 0.9, 0.95))
# meta.data_cols <- c("nGenes", "nUMIs", "MitoRatio", "Novelty")
# 
# # Loop through each cutoff and calculate cells passing QC
# for (i in 1:4){
#   
#   header_df <- data.frame(Sample = "Cutoffs", 
#                           gene_cutoff = as.numeric(cutoffs[[1]][i]),
#                           umi_cutoff = as.numeric(cutoffs[[2]][i]),
#                           mito_cutoff = as.numeric(cutoffs[[3]][i]),
#                           novelty_cutoff = as.numeric(cutoffs[[4]][i]))
#   
#   df <- data.frame(Sample = unique(raw_seurat@meta.data$Sample))
#   
#   for (j in 1:4){
#     
#     if(j == 3){  # for mito_cutoff, we need to keep cells below than cutoff
#       df1 <- raw_seurat@meta.data %>% 
#         dplyr::group_by(Sample) %>% 
#         dplyr::filter(!!rlang::sym(meta.data_cols[j]) <= cutoffs[[j]][i]) %>%
#         dplyr::count() %>%
#         dplyr::rename(Sample = identity(1), !!names(cutoffs)[j] :=  identity(2))
#     } else { # for all other cutoffs, we need to keep cells above than cutoff
#       df1 <- raw_seurat@meta.data %>% 
#         dplyr::group_by(Sample) %>% 
#         dplyr::filter(!!rlang::sym(meta.data_cols[j]) >= cutoffs[[j]][i]) %>%
#         dplyr::count() %>%
#         dplyr::rename(Sample = identity(1), !!names(cutoffs)[j] :=  identity(2))
#     }
#     df <- df %>% dplyr::left_join(df1, by=c("Sample"="Sample"))
#   }
#   summary_df <- dplyr::bind_rows(summary_df, header_df, df)
# }
# 
# # Save the summary
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb = wb, sheetName = paste0("Summary"))
# openxlsx::writeData(wb = wb, sheet = paste0("Summary"), x = summary_df)
# openxlsx::saveWorkbook(wb = wb, file = paste0(seurat_results, "Summary_QC.xlsx"), overwrite = TRUE)

# Choose appropriate cutoffs based on excel file we generated above. 
# NOTE: Increase Novelty score cutoff if there are cells with nUMIs > 50000
# Ambiguously demuxed cells are not assigned to any patient. So, remove cells 
# that dont map to any Patient.
gene_cutoff <- 250
umi_cutoff <- 500
mito_cutoff <- 0.2
ribo_cutoff <- 0.05
novelty_cutoff <- 0.8  	        # use 0.8 as starting point. Maximum 0.9
filtered_seurat <- base::subset(x = raw_seurat,
                                subset = (nGenes >= gene_cutoff) &
                                  # (RiboRatio >= ribo_cutoff) &
                                  # (!is.na(filtered_seurat@meta.data$Patient)) &
                                  (nUMIs >= umi_cutoff) & 
                                  (MitoRatio <= mito_cutoff) &
                                  (Novelty >= novelty_cutoff))

# Remove HTO assay from to avoid complications during integration, etc
if (demultiplexed == "yes"){
  filtered_seurat[["HTO"]] <- NULL
}

# Create whitelist file containing barcodes for each sample
if (whitelist == "yes"){
  
  # Extract barcodes and split by "_"
  bc <- filtered_seurat@meta.data$Cell
  
  barcodes <- data.frame(stringr::str_split_fixed(bc, "_", 2)) %>%
    dplyr::rename(Batch = identity(1), Barcodes = identity(2)) %>%
    dplyr::mutate(Barcodes = stringr::str_replace(Barcodes, "-1", ""))
  
  # Check how many barcodes are present in each batch
  barcodes %>% dplyr::group_by(Batch) %>% dplyr::count()
  
  # Save barcodes to individual csv files
  for (i in unique(barcodes$Batch)){
    whitelist <- barcodes %>% 
      dplyr::filter(Batch == i) %>% 
      dplyr::select(Barcodes)
    
    write.table(x = whitelist, 
                file = paste0(scripts_path, proj, "_", i, "_whitelist.csv"),
                row.names = FALSE,
                col.names = FALSE)
  }
  
  # # Read the file containing info to match SRR_ID to Batch
  # srr <- read.table(file = paste0(scripts_path, "SRR_Acc_Batch_Info.txt"))
  # colnames(srr) <- c("SRR_ID", "Batch")
  # 
  # # Save all barcodes in each batch as a separate csv file
  # for (i in unique(barcodes$Batch)){
  #   whitelist <- barcodes %>% filter(Batch == i) %>% select(Barcodes)
  #   # Match each batch with its SRR_Accession id
  #   prefix <- srr %>% filter(grepl(i, Batch))
  #   write.table(x = whitelist, 
  #               file = paste0(scripts_path, prefix$SRR_ID[1], "_", i, "_whitelist.csv"),
  #               row.names = FALSE,
  #               col.names = FALSE)
}

#****************************STEP 2C: SAVE THE DATA****************************#

# Create .rds object for filtered seurat object to load at any time
saveRDS(filtered_seurat, file=paste0(seurat_results, "filtered_seurat.rds"))

#*****************STEP 2D: VISUALIZE DATA BEFORE AND AFTER QC******************#

raw_metadata <- raw_seurat@meta.data
filtered_metadata <- filtered_seurat@meta.data

# Visualize the number of cell counts per sample
cell_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, fill=Sample)) + 
    geom_bar() +              
    theme_classic() +         #display with x and y axis lines and no gridlines
    labs(x = "Sample", y = "Number of Cells", title = stringr::str_wrap(paste0("Number of Cells ", tag), 30)) +
    #coord_cartesian(ylim = c(0, 20000)) +
    geom_text(stat="count", aes(label=after_stat(count)), vjust = -1) +
    my_theme
}

# Visualize the number of UMIs per cell
umi_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=nUMIs, fill=Sample)) +
    geom_violin() +        
    theme_classic() +       
    labs(x = "Sample", y = "Number of UMIs", title = stringr::str_wrap(paste0("Distribution of UMIs ", tag),30)) +
    coord_cartesian(ylim = c(100, 1000000)) +
    my_theme +        
    scale_y_log10(breaks = c(100, 1000, 10000, 100000, 1000000)) +  		     #display y axis in log scale
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = gene_cutoff, linetype = 2)
}

# Visualize the number of genes per cell
gene_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=nGenes, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "Number of Genes", title = stringr::str_wrap(paste0("Distribution of Genes ", tag),30)) +
    coord_cartesian(ylim = c(1, 30000)) +
    my_theme +        
    scale_y_log10() +    	
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = gene_cutoff, linetype = 2)
}

# Visualize the MitoRatio of each cell
mito_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=MitoRatio, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "MitoRatio", title = stringr::str_wrap(paste0("Distribution of MitoRatio ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = mito_cutoff, linetype = 2)
}

# Visualize the RiboRatio of each cell
ribo_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=RiboRatio, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "RiboRatio", title = stringr::str_wrap(paste0("Distribution of RiboRatio ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = ribo_cutoff, linetype = 2)
}

# Visualize the novelty or complexity of each cell
novelty_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=Novelty, fill = Sample)) +
    geom_violin() +     
    theme_classic() + 
    labs(x = "Sample", y = "Novelty Score", title = stringr::str_wrap(paste0("Distribution of Novelty Score ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = novelty_cutoff, linetype = 2)
}

# Visualize number of genes/cell, number of UMIs/cell & MitoRatio together.
# Bottom left quadrant : Poor quality cells with low genes & UMIs per cell 
# Top right quadrant   : Good quality cells with high genes & UMIs per cell
# Bottom right quadrant: Cells with low genes but high UMIs per cell. These 
# could be dying cells or population of low complexity cells (i.e erythrocytes)
gene_umi_mito_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=nUMIs, y=nGenes, color = MitoRatio)) +
    geom_point() +
    theme_classic() + 
    labs(x = "Number of UMIs", y = "Number of Genes",	 title = paste0("Distribution of UMIs, Genes & MitoRatio ", tag)) +
    my_theme + 
    coord_cartesian(xlim = c(100, 1000000), ylim = c(100, 20000)) +
    scale_x_log10(breaks = c(100, 1000, 10000, 100000, 1000000)) + 
    scale_y_log10() + 
    facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
    stat_smooth(method=lm, color="yellow") +
    geom_vline(xintercept = umi_cutoff) +    	#draw a vertical line at x=500 i.e.UMIs cutoff
    geom_hline(yintercept = gene_cutoff) +    #draw a horizontal line at y =250 i.e. Genes cutoff
    scale_color_viridis(option = "D", limits = c(0, 1)) 		# limits sets max and min values of gradient 
}

# Plot all QC metrics before and after QC
funcs <- c("cell_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "novelty_qc",
           "gene_umi_mito_qc")

filenames <- c("Cell_Counts", "UMI_Distribution", "Gene_Distribution",
               "MitoRatio_Distribution", "RiboRatio_Distribution", 
               "Novelty_Score_Distribution", "Genes_UMI_MitoRatio_Distribution")

for (i in 1:length(funcs)){
  
  # Plot QC metrics
  purrr::map2(.x = c("raw_metadata", "filtered_metadata"),
              .y = c("Pre QC", "Post QC"),
              .f = get(funcs[i])) %>% 
    cowplot::plot_grid(plotlist = .,
                       align = "hv",
                       axis = "tblr",
                       nrow = 2,  
                       ncol = 1, 
                       rel_widths = 1,
                       rel_heights = 1,
                       labels = NULL,
                       label_size = 14,
                       label_fontfamily = NULL,
                       label_fontface = "bold",
                       label_colour = NULL,
                       label_x = 0,
                       label_y = 1,
                       hjust = -0.5,
                       vjust = 1.5,
                       scale = 1,
                       greedy = TRUE,
                       byrow = TRUE,
                       cols = NULL,
                       rows = NULL)  
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("QC_", filenames[i], ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = seurat_results,
                  scale = 1,
                  width = dplyr::if_else(i==7, 17,8.5),
                  height = dplyr::if_else(i==7, 22, 11),
                  units = c("in"),	 
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
}

# ggplot2::ggplot(data = filtered_seurat@meta.data, aes(x=nUMIs, y=Sample, fill= after_stat(x))) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   viridis::scale_fill_viridis(name = "nUMIs", alpha = 1, begin = 0, end = 1, 
#                               direction = 1, discrete = FALSE, option = "D") +
#   labs(title = 'UMI Distribution') +
#   coord_cartesian(xlim = c(100, 10000)) +
#   scale_x_continuous(breaks = c(100, 1000, 2500, 5000, 10000))+
#   #hrbrthemes::theme_ipsum() +
#   theme(legend.position="right",
#         panel.spacing = unit(0.1, "lines"),
#         strip.text.x = element_text(size = 8))
# 
# ggsave("1.jpg", bg="white")  

#******************************************************************************#
#                       STEP 3: RUN THE STANDARD PIPELINE                      #
#******************************************************************************#

# Load rds file of filtered_seurat object
filtered_seurat <- base::readRDS(paste0(seurat_results, "filtered_seurat.rds"))

celltype <- NULL
res <- "integrated_snn_res.1.4"
filt_data <- filtered_seurat
sct_data <- sctransform_data(filt_data)
plot_pre_integration(sct_data)
integ_data <- integrate_data(sct_data, filt_data)
integ_data <- cluster_data(integ_data)
integ_data <- clean_data(integ_data, res)
save_data(integ_data)
plot_post_integration(integ_data, res)
plot_modules(res, "All Markers")

#******************************************************************************#
