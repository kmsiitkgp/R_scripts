# Visit this to see video of all steps in GSEA
# https://www.youtube.com/watch?v=Bzu4_yDcBLY

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
library("TCGAbiolinks")         
library("ensembldb")            
library("AnnotationHub")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("clusterProfiler")
library("fgsea")
library("enrichplot")
library("pathview")
library("msigdbr")

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")
library("stringr")

# Graph plotting packages
library("ggplot2")
library("viridis")
library("RColorBrewer")
library("ggrepel")
library("ggbeeswarm")

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

parent_path <- "C:/Users/KailasammS/Desktop/"
results_path <- "C:/Users/KailasammS/Desktop/"

# Check MSig database version
packageVersion("msigdbr")

# Choose which gene-set category and subcategory to use for analysis by going 
# through msig_collections dataframe. Most common ones are C5_GO:BP and H.
msigdb_collections <- msigdbr::msigdbr_collections()
msigdb_species <- msigdbr::msigdbr_species()
msigdb_genesets <- list(cat = c("H", "C5", "C5"),
                        subcat = c("", "GO:BP", "GO:MF"))

# Choose species
species <- "Homo sapiens"
species <- "Mus musculus"

# Choose which database you want to get gene annotations from.
# NOTE: Most RNA Seq data has Ensembl ids as rownames
database <- "Entrez"
database <- "Ensembl"

# Choose database for encrichGO()
species_db <- dplyr::if_else(species == "Mus musculus", "org.Mm.eg.db", "org.Hs.eg.db")

# Choose whether to plot GSEA plot
plot_fGSEA <- FALSE

#******************************************************************************#
#                            STEP 1: GET ANNOTATIONS                           #
#******************************************************************************#

get_annotations <- function(database){
  # NOTE: Before running FindAllMarkers(), get gene description from ensembl
  
  # To convert from Ensembl to gene symbol, you can use AnnotationHub/AnnotationDbi
  # AnnotationHub seems to be more comprehensive for Ensembl-gene symbol conversion.
  # However, to convert from Entrez to gene symbol, you can ONLY use AnnotationDbi
  
  #**************************GET ENSEMBL ANNOTATIONS***************************#
  # Connect to AnnotationHub
  ah <- AnnotationHub()
  
  # Access the Ensembl database for organism
  ahDb <- AnnotationHub::query(ah,
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
    dplyr::select(gene_id, entrezid, gene_name, seq_name, gene_seq_start, gene_seq_end, gene_biotype, description) %>%
    dplyr::rename(ensembl_id = gene_id, entrez_id = entrezid, chr = seq_name, start = gene_seq_start, end = gene_seq_end) %>%
    dplyr::distinct(ensembl_id, .keep_all = TRUE)
  
  #***************************GET ENTREZ ANNOTATIONS***************************# 
  
  # # Map Entrez id with "ENSEMBL" or "SYMBOL" whichever maps highest genes
  # # Here we map with SYMBOL
  # mapping <- AnnotationDbi::mapIds(x = org.Hs.eg.db, keys = read_data$ENTREZID,
  #                                  keytype = "ENTREZID", column = "SYMBOL")
  # 
  # mapping <- as.data.frame(do.call(cbind, list(mapping))) %>%
  #   tibble::rownames_to_column("ENTREZID")
  # 
  # colnames(mapping) <- c("ENTREZID","GENE")
  # 
  # # Merge read_data with mapped data
  # read_data <- read_data %>%
  #   dplyr::left_join(mapping, by=c("ENTREZID"="ENTREZID")) %>%
  #   dplyr::select(GENE, everything(), -ENTREZID)
  # 
  # # Check how many ENTREZIDs couldn't be mapped
  # unmapped <- nrow(read_data %>% dplyr::filter(is.na(.[[1]])))
  # cat("Genes in total :", nrow(read_data), "\t")
  # cat("Genes mapped   :", nrow(read_data)-unmapped, "\t")
  # cat("Genes unmapped :", unmapped, "\t")
  
  if (species == "Homo sapiens"){
    entrez <- AnnotationDbi::select(x = org.Hs.eg.db, keys = keys(org.Hs.eg.db),
                                    keytype = "ENTREZID", column = "SYMBOL")
  } else{
    entrez <- AnnotationDbi::select(x = org.Mm.eg.db, keys = keys(org.Mm.eg.db),
                                    keytype = "ENTREZID", column = "SYMBOL")
  }
  colnames(entrez) <- c("ID", "SYMBOL")
  
  # Choose which annotations to use
  if(database == "Ensembl"){
    annotations <- ensembl
  } else{
    annotations <- entrez
  }
  
  return(annotations)
  DBI::dbDisconnect(conn = ah)
}

annotations <- get_annotations(database)

#******************************************************************************#
#        STEP 2: GET REQUIRED GENE SET COLLECTIONS FROM BROAD DATABASE         #
#******************************************************************************#

for (i in 1:length(msigdb_genesets$cat)){
  msigdb_collection <- msigdbr::msigdbr(species = species,
                                        category = msigdb_genesets$cat[i],        
                                        subcategory = msigdb_genesets$subcat[i]) 
  
  msigdb_collection <- msigdb_collection %>% 
    dplyr::select(gs_name, ensembl_gene)
  
  # Assign the object to its corresponding variable
  assign(x = paste0(msigdb_genesets$cat[i], "_", msigdb_genesets$subcat[i]),
         value = msigdb_collection)
}

#******************************************************************************#  
#             OVER-REPRESENTATION ANALYSIS USING CLUSTER PROFILER              #
#******************************************************************************#

ora <- function(DEGs_df){
  
  # Over-representation analysis looks at whether a list of genes that you have
  # already separated out (i.e. identified as DEGS) associate significantly with
  # certain pathways. It doesn't consider the expression level of the genes.
  
  # We use ONLY SIGNIFICANT DEGS for this analysis. 
  # Filter out genes with padj > 0.05 and abs(log2FoldChange) >= 0.58
  DEGs_df_sig <- DEGs_df %>% 
    dplyr::filter(padj < 0.05 & !is.na(padj)) %>%
    dplyr::left_join(annotations, by=c("SYMBOL"="gene_name")) %>%
    dplyr::distinct_at("ensembl_id", .keep_all = TRUE)
  
  # Run enricher() on Hallmark gene sets
  enriched_hallmark <- clusterProfiler::enricher(gene = DEGs_df_sig$ensembl_id,
                                                 pvalueCutoff = 0.05,
                                                 pAdjustMethod = "BH",
                                                 universe = NULL,
                                                 minGSSize = 10,
                                                 maxGSSize = 500,
                                                 qvalueCutoff = 0.2,
                                                 TERM2GENE = H_, #collection
                                                 TERM2NAME = NA)
  
  # Run enrichGO() on GO BP gene sets
  enriched_GO_BP <- clusterProfiler::enrichGO(gene = DEGs_df_sig$entrez_id,
                                              OrgDb = species_db,     
                                              keyType = "ENTREZID",   #id
                                              ont = "BP",  #"MF" or "CC"
                                              pvalueCutoff = 0.05,
                                              pAdjustMethod = "BH",
                                              universe = NULL,
                                              minGSSize = 10,
                                              maxGSSize = 500,
                                              readable = FALSE,
                                              pool = FALSE)
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  
  for (type in c("enriched_hallmark", "enriched_GO_BP")){
    
    # Format results to make it plottable
    enriched_result <- get(type)@result %>% 
      tibble::remove_rownames() %>%
      dplyr::select(everything(), -c("ID")) %>%
      tidyr::separate(col = "GeneRatio", into = c("DEGs_df.pathway", "DEGs_df.collection")) %>%
      tidyr::separate(col = "BgRatio", into = c("Genes.pathway", "Genes.collection")) %>%
      dplyr::mutate_at(c("DEGs_df.pathway", "DEGs_df.collection", "Genes.pathway", "Genes.collection"), as.numeric) %>%
      dplyr::mutate(k.K = DEGs_df.pathway/Genes.pathway) %>%
      dplyr::filter(p.adjust < 0.05) %>%
      dplyr::arrange(desc(k.K)) %>% 
      dplyr::mutate(Description = gsub("HALLMARK_", "", Description),
                    Description = gsub("_", " ", Description),
                    # Description = gsub("ENDOPLASMIC RETICULUM", "ER", Description),
                    # Description = gsub("-Positive", "+", Description),
                    # Description = gsub("Nadh", "NADH", Description),
                    #length = stringr::str_length(Description),
                    Description = stringr::str_trunc(Description, 45, "right"),
                    Description = stringr::str_to_upper(Description),
                    Description = stringr::str_wrap(Description, width = 22))
    
    # Visualize results if there are significant pathways. Else, skip plotting.
    # It is difficult to control the width of bars in ggplot. Since, we plot
    # top 12 pathways, we insert dummy entries to make the data frame have 12
    # pathways "if" the data frame has less than 12 pathways
    
    if(nrow(enriched_result) > 0){
      
      if(nrow(enriched_result) < 12){
        nrows <- nrow(enriched_result)
        enriched_result[(nrows+1):12,] <- seq(nrows+1:(12-nrows))
        enriched_result$p.adjust[(nrows+1):12]  <- rep(c(0), each = 12-nrows)
        enriched_result$k.K[(nrows+1):12]  <- rep(c(0), each = 12-nrows)
      } else{
        enriched_result <- enriched_result %>% 
          dplyr::slice_max(k.K, n = 12)
      }
      
      # Read the corrected file after adjusting TNFa, KRAS etc
      #enriched_result <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "GSEA_plots/(C5) GSEA results.xlsx"))
      
      # Plot bar plots
      # (NOT RECOMMEDNED since it gives info only on enrichment ratio & pvalue)
      ggplot2::ggplot(data = enriched_result,
                      aes(x = k.K,
                          y = Description,
                          fill = p.adjust)) +
        ggplot2::geom_col(width = 0.75) +
        ggplot2::theme_classic() +
        ggplot2::labs(x = "Enrichment Ratio",
                      y = "",
                      fill = "padj",
                      title = "Pathways") +
        #ggplot2::coord_cartesian(xlim = c(0, max(abs(enriched_result$k.K)))) +
        ggplot2::theme(aspect.ratio = 2,
                       plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
                       axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
                       legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0, vjust = 1),
                       legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                       legend.position = "bottom", #"right",
                       legend.justification = "center",
                       legend.direction = "horizontal", #"vertical",
                       legend.key.height = unit(0.5, 'cm'),     #unit(0.75, 'cm'),
                       legend.key.width = unit(1.25, 'cm')) +   #unit(0.5, 'cm'),
        viridis::scale_fill_viridis(option = "viridis", limits = c(0, 0.05))
      
      # Save the plot
      ggplot2::ggsave(filename = paste0("(ORA)_Bar_plot", type, "_.tiff"),
                      plot = last_plot(),
                      device = "jpeg",
                      path = results_path,
                      scale = 1,
                      width = 6,
                      height = 7,
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = NULL)
      
      # Plot dot plots (RECOMMENDED)
      ggplot2::ggplot(data = enriched_result,
                      aes(x = k.K,
                          y = Description,
                          #label = Upstream.Regulator,
                          size = Count,
                          color = p.adjust)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Enrichment Ratio",
                      y = "",
                      title = "",
                      fill = "padj") +
        ggplot2::theme(plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
                       axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
                       legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0, vjust = 1),
                       legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                       #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                       legend.position = "right",
                       legend.justification = "left",
                       legend.direction = "vertical",
                       legend.key.height = unit(0.5, 'cm'),
                       legend.key.width = unit(0.5, 'cm')) +
        viridis::scale_color_viridis(option = "viridis") +
        ggplot2::scale_size(breaks = sapply(as.vector(quantile(c(min(enriched_result$Count), max(enriched_result$Count)))), floor))
      
      # Save the plot
      ggplot2::ggsave(filename = paste0("(ORA)_Dot_plot", type, "_.tiff"),
                      plot = last_plot(),
                      device = "jpeg",
                      path = results_path,
                      scale = 1,
                      width = 6,
                      height = 7,
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = NULL)
      
      openxlsx::addWorksheet(wb, sheetName = paste0("(ORA)_", type))
      openxlsx::writeData(wb, sheet = paste0("(ORA)_", type), x = enriched_result, 
                          rowNames = FALSE)
    } else { print("No significant gene sets")}
  }
  
  # Save the results in excel file
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "_ORA_results.xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#                          GSEA ANALYSIS USING FGSEA                           #
#******************************************************************************#

# Enrichment analysis takes differential data from every measured gene and looks
# for pathways displaying significantly coordinated shifts in those values.

fgsea <- function(DEGs_df, collection){
  
  #***************************DEFINE stats PARAMETER***************************#
  
  # NOTE: ALL genes MUST be used for this analysis, NOT just DEGs. 
  # NOTE: Genes MUST be sorted in descending fold change.
  # NOTE: Genes MUST be stored in list format, not as a dataframe.
  colnames(DEGs_df)[1] <- "Ensembl"
  DEGs_df <- DEGs_df %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(desc(log2FoldChange))
  
  DEGs_list <- DEGs_df$log2FoldChange
  names(DEGs_list) <- DEGs_df$Ensembl
  
  #*************************DEFINE scoreType PARAMETER*************************#
  
  # Define score type in fgseaMultilevel() based on fold change values.
  # NOTE: Use "pos", if you are ONLY interested in activated pathways.
  # NOTE: Use "neg", if you are ONLY interested in inhibited pathways. 
  # NOTE: Else, use "std" for both activated & inhibited pathways. 
  score_type <- dplyr::if_else(max(DEGs_list) > 0 & min(DEGs_list) < 0, "std", 
                               dplyr::if_else(max(DEGs_list) < 0 & min(DEGs_list) < 0, "neg", "pos"))
  
  #*************************DEFINE pathways PARAMETER**************************#
  
  # NOTE: Unlike clusterProfiler::GSEA(), fgsea::fgseaMultilevel() needs gene 
  # sets in a specific list format
  fgsea_pathways <- list()
  
  for (pathway in unique(get(collection)$gs_name)){
    l1 <- as.list(get(collection) %>% 
                    dplyr::filter(gs_name == pathway) %>% 
                    dplyr::select(ensembl_gene) %>%
                    dplyr::distinct())
    
    names(l1) <- pathway 
    fgsea_pathways <- append(fgsea_pathways, l1)
  }
  
  #**********************************RUN fGSEA*********************************#
  
  fgsea_results <- fgsea::fgseaMultilevel(pathways = fgsea_pathways,
                                          stats = DEGs_list,
                                          scoreType = score_type,
                                          sampleSize = 101,
                                          minSize = 1,
                                          maxSize = length(DEGs_list) - 1,
                                          eps = 1e-50,
                                          nproc = 0,
                                          gseaParam = 1,
                                          BPPARAM = NULL,
                                          nPermSimple = 1000)
  
  #*******************************FORMAT RESULTS*******************************#
  
  # NOTE: Output of fgsea is a data.table & data.frame). "leadingEdge" column
  # is a list of genes. So, DO NOT FORCE the output of fgsea to a dataframe as 
  # this will lead to data loss from "leadingEdge" column & affect plotting 
  # using fgsea::plotEnrichment()
  
  # Filter out non-significant pathways
  fgsea_results <- fgsea_results %>% 
    dplyr::filter(padj < 0.05) %>%
    dplyr::mutate(abs_NES = abs(NES)) %>%
    dplyr::arrange(desc(abs_NES))
  
  # Identify overlapping pathways and collpase them into major pathways
  concise_fgsea <- fgsea::collapsePathways(fgseaRes = fgsea_results,
                                           pathways = fgsea_pathways,
                                           stats = DEGs_list)
  # Filter out overlapping pathways
  fgsea_results <- fgsea_results %>% 
    dplyr::filter(pathway %in% concise_fgsea$mainPathways) %>%
    dplyr::mutate(direction = dplyr::if_else(NES > 0, "Activated", "Inhibited"))
  
  # Convert ensembl ids to gene symbol
  max_len <- max(unlist(lapply(X=fgsea_results$leadingEdge, FUN=length)))
  genes_df <- data.frame(matrix(NA, nrow=max_len))
  
  for (i in 1:nrow(fgsea_results)){
    l1 <- annotations %>% 
      dplyr::filter(ensembl_id %in% fgsea_results$leadingEdge[[i]]) %>% 
      dplyr::select(gene_name) %>% 
      unlist(., use.names=FALSE)
    l1 <- c(l1, rep(x=NA, times=max_len-length(l1)))
    
    genes_df <- dplyr::bind_cols(genes_df, l1)
    colnames(genes_df)[i+1] <- fgsea_results$pathway[i]
  }
  genes_df <- genes_df[,-1]
  
  # Save the results in excel file
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = gsub(":", "_", collection))
  openxlsx::writeData(wb, sheet = gsub(":", "_", collection), x = fgsea_results, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Genes")
  openxlsx::writeData(wb, sheet = "Genes", x = genes_df, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "fGSEA_Results.xlsx"),
                         overwrite = TRUE)
  
  #*******************PLOT ENRICHMENT PLOT FOR EACH PATHWAY********************#
  
  if (plot_fGSEA){
    for (p in concise_fgsea$mainPathways){
      fgsea::plotEnrichment(pathway = fgsea_results %>%
                              dplyr::filter(pathway == p) %>%
                              dplyr::select(leadingEdge) %>%
                              unlist(., use.names = FALSE),
                            stats = DEGs_list,
                            gseaParam = 1,
                            ticksSize = 0.2)
      
      ggplot2::ggsave(filename = paste0("(fGSEA)_Enrichment_plot", gsub(":", "_", collection),"_", p, ".tiff"),
                      plot = last_plot(),
                      device = "jpeg",
                      path = results_path,
                      scale = 1,
                      width = 6,
                      height = 7,
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = NULL)
    }
  }
  
  #****************PLOT SUMMARY OF ENRICHED PATHWAYS AS DOT PLOT***************#
  
  # Visualize results if there are significant pathways. Else, skip plotting
  if(nrow(fgsea_results) > 0){
    
    # It is difficult to control the width of bars in ggplot. Since, we plot
    # top 12 pathways, we insert dummy entries to make the data frame have 12
    # pathways "if" the data frame has less than 12 pathways
    if (nrow(fgsea_results) < 12){
      nrows <- nrow(fgsea_results)
      fgsea_results[(nrows+1):12,] <- seq(nrows+1:(12-nrows))
      fgsea_results$NES[(nrows+1):12]  <- rep(c(0), each = 12-nrows)
      fgsea_results$direction[(nrows+1):12] <- rep(c("Inhibited"), each = 12-nrows)
    } else {
      fgsea_results <- fgsea_results %>% 
        dplyr::slice_max(abs_NES, n = 12)
    }
    
    # Modify pathway names to make the plot pretty
    fgsea_results_pretty <- fgsea_results %>%
      dplyr::mutate(pathway = gsub("HALLMARK_|SA_|SIG_|NABA_|GOBP_|GOMF_", "", pathway),
                    pathway = gsub("_", " ", pathway),
                    pathway = gsub("ENDOPLASMIC RETICULUM", "ER", pathway),
                    #pathway = stringr::str_trunc(pathway, 45, "right"),
                    #pathway = stringr::str_to_title(pathway),
                    #length = stringr::str_length(pathway),
                    pathway = stringr::str_wrap(pathway, width = 22))
    
    ggplot2::ggplot(data = fgsea_results_pretty, 
                    aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
      # fill = direction means direction will be arranged in alphabetical order.
      # So, if you had labeled direction as "Upregulated" and "Downregulated",
      # then first color in scale_fill_manual() will be assigned to 
      # "downregulated" and it will be labeled as "Activated in Males" in the
      # plot. So, be careful.
      ggplot2::geom_col(width = 0.75) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Normalized Enrichment Score(NES)",
                    y = "",
                    title = "GSEA",
                    fill = "") +
      ggplot2::coord_cartesian(xlim = c(floor(-max(abs(fgsea_results$NES), na.rm=TRUE)), ceiling(max(abs(fgsea_results$NES), na.rm=TRUE)))) +
      ggplot2::theme(#aspect.ratio = 2,
        plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
        axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
        legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
        legend.position = "bottom",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.25, 'cm')) +
      ggplot2::scale_fill_manual(labels=c("Activated", "Inhibited"),
                                 values = c(RColorBrewer::brewer.pal(11, "RdYlBu")[c(1)],
                                            RColorBrewer::brewer.pal(11, "RdYlBu")[c(11)]))
    
    # Save the plot
    ggplot2::ggsave(filename = paste0("(fGSEA)_Dot_plot", gsub(":", "_", collection), ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = results_path,
                    scale = 1,
                    width = 6,
                    height = 7,
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = NULL)
  } else { print("No significant gene sets")}
}

#******************************************************************************#
#                     GSEA ANALYSIS USING CLUSTER PROFILER                     #
#******************************************************************************#

# Enrichment analysis takes differential data from every measured gene and looks
# for pathways displaying significantly coordinated shifts in those values.

gsea <- function(DEGs_df, collection){  
  
  #*************************DEFINE geneList PARAMETER**************************#
  
  # NOTE: ALL genes MUST be used for this analysis, NOT just DEGs. 
  # NOTE: Genes MUST be sorted in descending fold change.
  # NOTE: Genes MUST be stored in list format, not as a dataframe.
  colnames(DEGs_df)[1] <- "Ensembl"
  DEGs_df <- DEGs_df %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(desc(log2FoldChange))
  
  DEGs_list <- DEGs_df$log2FoldChange
  names(DEGs_list) <- DEGs_df$Ensembl
  
  #**********************************RUN GSEA**********************************#
  
  gsea_results <- clusterProfiler::GSEA(geneList = DEGs_list,
                                        exponent = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        eps = 1e-10,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH",
                                        TERM2GENE = get(collection),
                                        TERM2NAME = NA,
                                        verbose = TRUE,
                                        seed = FALSE,
                                        by = "fgsea")
  
  #*******************************FORMAT RESULTS*******************************#
  
  # Filter out non-significant pathways
  gsea_results <- gsea_results@result %>% 
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::mutate(abs_NES = abs(NES)) %>%
    dplyr::arrange(desc(abs_NES))
  
  # Filter out overlapping pathways
  gsea_results <- gsea_results %>% 
    dplyr::mutate(direction = dplyr::if_else(NES > 0, "Activated", "Inhibited"))
  
  # Convert ensembl ids to gene symbol
  max_len <- max(unlist(lapply(X=stringr::str_split(string = gsea_results$core_enrichment, pattern = "/"), FUN=length)))
  genes_df <- data.frame(matrix(NA, nrow=max_len))
  
  for (i in 1:nrow(gsea_results)){
    l1 <- annotations %>% 
      dplyr::filter(ensembl_id %in% unlist(stringr::str_split(string = gsea_results$core_enrichment[[i]], pattern = "/"))) %>% 
      dplyr::select(gene_name) %>% 
      unlist(., use.names=FALSE)
    l1 <- c(l1, rep(x=NA, times=max_len-length(l1)))
    
    genes_df <- dplyr::bind_cols(genes_df, l1)
    colnames(genes_df)[i+1] <- gsea_results$Description[i]
  }
  genes_df <- genes_df[,-1]
  
  # Save the results in excel file
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = gsub(":", "_", collection))
  openxlsx::writeData(wb, sheet = gsub(":", "_", collection), x = gsea_results, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Genes")
  openxlsx::writeData(wb, sheet = "Genes", x = genes_df, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "GSEA_Results.xlsx"),
                         overwrite = TRUE)
  
  #*******************PLOT ENRICHMENT PLOT FOR EACH PATHWAY********************#
  
  #****************PLOT SUMMARY OF ENRICHED PATHWAYS AS DOT PLOT***************#
  
  # Visualize results if there are significant pathways. Else, skip plotting
  if(nrow(gsea_results) > 0){
    
    # It is difficult to control the width of bars in ggplot. Since, we plot
    # top 12 pathways, we insert dummy entries to make the data frame have 12
    # pathways "if" the data frame has less than 12 pathways
    if (nrow(gsea_results) < 12){
      nrows <- nrow(gsea_results)
      gsea_results[(nrows+1):12,] <- seq(nrows+1:(12-nrows))
      gsea_results$NES[(nrows+1):12]  <- rep(c(0), each = 12-nrows)
      gsea_results$direction[(nrows+1):12] <- rep(c("Inhibited"), each = 12-nrows)
    } else {
      gsea_results <- gsea_results %>% 
        dplyr::slice_max(abs_NES, n = 12)
    }
    
    # Modify pathway names to make the plot pretty
    gsea_results_pretty <- gsea_results %>%
      dplyr::mutate(pathway = gsub("HALLMARK_|SA_|SIG_|NABA_|GOBP_|GOMF_", "", Description),
                    pathway = gsub("_", " ", pathway),
                    pathway = gsub("ENDOPLASMIC RETICULUM", "ER", pathway),
                    #pathway = stringr::str_trunc(pathway, 45, "right"),
                    #pathway = stringr::str_to_title(pathway),
                    #length = stringr::str_length(pathway),
                    pathway = stringr::str_wrap(pathway, width = 22))
    
    ggplot2::ggplot(data = gsea_results_pretty, 
                    aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
      # fill = direction means direction will be arranged in alphabetical order.
      # So, if you had labeled direction as "Upregulated" and "Downregulated",
      # then first color in scale_fill_manual() will be assigned to 
      # "downregulated" and it will be labeled as "Activated in Males" in the
      # plot. So, be careful.
      ggplot2::geom_col(width = 0.75) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Normalized Enrichment Score(NES)",
                    y = "",
                    title = "GSEA",
                    fill = "") +
      ggplot2::coord_cartesian(xlim = c(floor(-max(abs(gsea_results$NES), na.rm=TRUE)), ceiling(max(abs(gsea_results$NES), na.rm=TRUE)))) +
      ggplot2::theme(#aspect.ratio = 2,
        plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
        axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
        legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
        legend.position = "bottom",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.25, 'cm')) +
      ggplot2::scale_fill_manual(labels=c("Activated", "Inhibited"),
                                 values = c(RColorBrewer::brewer.pal(11, "RdYlBu")[c(1)],
                                            RColorBrewer::brewer.pal(11, "RdYlBu")[c(11)]))
    
    # Save the plot
    ggplot2::ggsave(filename = paste0("(GSEA)_Dot_plot", gsub(":", "_", collection), ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = results_path,
                    scale = 1,
                    width = 6,
                    height = 7,
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = NULL)
  } else { print("No significant gene sets")}
}

#******************************************************************************#
#                           STEP 3: PERFORM ANALYSIS                           #
#******************************************************************************#

# The data frame we use as input to enricher() MUST have ONLY 2 columns: 
# one containing the genes and other containing the corresponding gene sets.
# The list of genes we specify in parameter "gene" within enricher() will be 
# compared against the genes in the input dataframe. So, make sure they are
# in the same format (gene name/ensembl_id/etc)

# Read DEGs_df
# NOTE: SYMBOL, padj, log2FoldChange columns MUST be present
data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Results__id_Y_pos_vs_Y_neg_DEGs.xlsx"))
colnames(data)[1] <- "Ensembl" 

g <- g %>% dplyr::select("Intersection#7", "Intersection#10")
ge <- c(g[,1], g[,2])
ge <- ge[!is.na(ge)]
data <- data %>% dplyr::filter(SYMBOL %in% ge)
collection <- c("C5_GO:BP")

ora(data)  
fgsea(data, collection)
gsea(data, collection)

#******************************************************************************#
#         (OPTIONAL): PLOT IPA CANONICAL PATHWAYS & UPSTREAM REGULATORS        #
#******************************************************************************#

parent_path <- "C:/Users/KailasammS/Desktop/"

# Re-plot Canonical Pathways from IPA
for (celltype in c("Epithelial", "Myeloid", "Fibroblasts", "Lymphoid")){
  
  ipa_pathways <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Canonical_", celltype, ".xlsx"),
                                      startRow = 2)
  
  colnames(ipa_pathways) <- c("Pathways", "Neg.log10p", "Ratio", "zscore", "Molecules")
  
  ipa_pathways <- ipa_pathways %>% 
    dplyr::filter(Neg.log10p > 1.301) %>%
    dplyr::mutate(direction = dplyr::if_else(zscore > 0, "Activated in Males", "Inhibited in Males"),
                  Pathways = gsub("Reactive Oxygen Species", "ROS", Pathways),
                  Pathways = stringr::str_trunc(Pathways, 50, "right"),
                  Pathways = stringr::str_wrap(Pathways, width = 30),
                  abs_zscore = abs(zscore)) %>%
    dplyr::arrange(desc(abs_zscore))
  
  if (celltype == "Epithelial"){
    ipa_pathways <- ipa_pathways %>% dplyr::slice(1,2,4,5,6,7,9,10,11)
  }
  if (celltype == "Myeloid"){
    ipa_pathways <- ipa_pathways %>% dplyr::slice(1,3,4,10,11,12,13,15,16,18)
  }
  if (celltype == "Fibroblasts"){
    ipa_pathways <- ipa_pathways %>% dplyr::slice(1,4,6,8,11,12,13,15,16,18)
  }
  if (celltype == "Lymphoid"){
    ipa_pathways <- ipa_pathways %>% dplyr::slice(1,2,5,6,7,9,10,11,12,14)
  }
  
  # Visualize results if there are significant pathways. Else, skip plotting
  if(nrow(ipa_pathways) > 0){
    
    # Determine limits to be used in plots
    max_x <- max(ipa_pathways$Neg.log10p)
    
    # It is difficult to control the width of bars in ggplot. Since, we plot
    # top 10 pathways, we insert dummy entries to make the data frame have 10
    # pathways "if" the data frame has less than 10 pathways
    if (nrow(ipa_pathways) < 10){
      nrows <- nrow(ipa_pathways)
      ipa_pathways[(nrows+1):10,] <- seq(nrows+1:(10-nrows))
      ipa_pathways$direction[(nrows+1):10]  <- rep(c("Activated in Males"), each = 10-nrows)
      ipa_pathways$abs_zscore[(nrows+1):10]  <- rep(c(0), each = 10-nrows)
      ipa_pathways$Neg.log10p[(nrows+1):10]  <- rep(c(100), each = 10-nrows)
    } else {
      ipa_pathways <- ipa_pathways %>% 
        dplyr::slice_max(abs_zscore, n = 10)
    }
    
    # Plot column graph
    ggplot2::ggplot(data = ipa_pathways, 
                    aes(x = Neg.log10p, y = reorder(Pathways, Neg.log10p), fill = direction)) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::scale_x_continuous(limits = c(0, ceiling(max_x)),
                                  breaks = seq(from=0, to=ceiling(max_x), length.out=5),
                                  position = "bottom",                   # X-axis plot label location
                                  expand = c(0, 0)) +                    # Remove padding on left and right   
      ggplot2::scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) + # Give 0.5 padding on top & bottom
      ggplot2::labs(x = " -log10(p-value)",
                    y = "",
                    title = paste0(""),
                    fill = "") +
      ggplot2::theme_classic() +
      ggplot2::theme(#panel.background = element_rect(fill = "white"), 
        panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3), 
        axis.ticks.length = unit(0, "mm"),
        plot.title =   element_text(family="sans", face="plain", colour="black", size=15, hjust = 0.5),
        plot.caption = element_text(family="sans", face="plain", colour="black", size=15, hjust = 0),
        axis.title.x = element_text(family="sans", face="plain", colour="black", size=15, hjust = 0.5),
        axis.title.y = element_text(family="sans", face="plain", colour="black", size=15, hjust = 0.5),
        axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
        legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
        #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
        legend.position = "bottom",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm')) +
      ggplot2::scale_fill_manual(labels=c("Activated in Males", "Inhibited in Males"),
                                 values = c(RColorBrewer::brewer.pal(11, "RdYlBu")[c(1)],
                                            RColorBrewer::brewer.pal(11, "RdYlBu")[c(11)]))
    # Save the plot
    ggplot2::ggsave(filename = paste0("(IPA_CP) ", celltype, ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = parent_path,
                    scale = 1,
                    #width = 7.5,
                    #height = 6,
                    units = c("in"),
                    dpi = 600,
                    limitsize = FALSE,
                    bg = NULL)
  }
}

# Re-plot Upstream Regulators from IPA
for (celltype in c("Epithelial", "Myeloid", "Fibroblasts", "Lymphoid")){
  
  upstream <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Upstream_", celltype, ".xlsx"),
                                  startRow = 2)
  
  # Remove - from colnames as R considers them as minus when using dplyr
  colnames(upstream) <- gsub("-", "", colnames(upstream))
  
  # Filter out unwanted molecules & convert Expr.Log.Ratio & Activation.zscore to 
  # numbers as they are stored as character
  upstream_subset <- upstream %>% 
    dplyr::mutate_at(c("Expr.Log.Ratio", "Activation.zscore"), as.numeric) %>% 
    dplyr::filter(abs(Activation.zscore) >= 2 & 
                    pvalue.of.overlap < 0.05 & 
                    abs(Expr.Log.Ratio) >= 0.58 &
                    !is.na(Expr.Log.Ratio)) %>%
    dplyr::mutate(Molecule.Type = gsub("ligand-dependent nuclear |transcription |translation |transmembrane |G-protein coupled ", "", Molecule.Type),
                  Molecule.Type = gsub("peptidase|phosphatase|kinase|enzyme", "Enzymes", Molecule.Type),
                  Molecule.Type = gsub("ion channel|transporter|receptor", "Receptors & Transporters", Molecule.Type),
                  Molecule.Type = gsub("cytokine|growth factor", "Cytokines & Growth factors", Molecule.Type),
                  Molecule.Type = gsub("regulator", "Transcription & Translation Regulators", Molecule.Type),
                  Molecule.Type = gsub("other|complex|group", "Other Regulators",  Molecule.Type))
  
  # Determine limits to be used in plots
  max_y <- max(abs(upstream_subset$Activation.zscore))
  min_y <- -max(abs(upstream_subset$Activation.zscore))
  max_x <- max(abs(upstream_subset$Expr.Log.Ratio))
  min_x <- -max(abs(upstream_subset$Expr.Log.Ratio))
  
  # Plot scatter plot
  ggplot2::ggplot(data = upstream_subset, 
                  mapping = aes(x = Expr.Log.Ratio, 
                                y = Activation.zscore, 
                                color = Predicted.Activation.State,
                                label = Upstream.Regulator)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(size = 3.5,
                             # To give space between points
                             #point.padding = 0.5,
                             # To give space between labels
                             #label.padding = 0.5,
                             box.padding = 0.5,
                             # To ALWAYS show all labels, even overlapping labels
                             max.overlaps = Inf,
                             # Repel away from the left & right edges.
                             xlim = c(NA, NA),
                             # Do not repel from bottom edges, but repel from top edge.
                             ylim = c(-Inf,NA),
                             # To draw any line segments. Set to 0 to ALWAYS draw
                             min.segment.length = 1000, #0.5,
                             # VERY IMPORTANT for proper positioning of labels
                             position = position_quasirandom()) +
    ggplot2::labs(x = "log2 Fold Change",
                  y = "Activation z-score",
                  title = "",
                  fill = "Direction") +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 2,
                   plot.title =   element_text(family="sans", face="bold", colour="black", size=18, hjust = 0.5),
                   plot.caption = element_text(family="sans", face="bold", colour="black", size=18, hjust = 0),
                   axis.title.x = element_text(family="sans", face="bold", colour="black", size=18, hjust = 0.5),
                   axis.title.y = element_text(family="sans", face="bold", colour="black", size=18, hjust = 0.5),
                   axis.text.x =  element_text(family="sans", face="bold", colour="black", size=18, hjust = 0.5),
                   axis.text.y =  element_text(family="sans", face="bold", colour="black", size=18, hjust = 1),
                   legend.title = element_text(family="sans", face="bold", colour="black", size=18, hjust = 0.5),
                   legend.text =  element_text(family="sans", face="bold", colour="black", size=18, hjust = 0.5),
                   #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                   legend.position = "bottom",
                   legend.justification = "left",
                   legend.direction = "horizontal",
                   legend.key.height = unit(0.5, 'cm'),
                   legend.key.width = unit(0.5, 'cm')) +
    ggplot2::coord_cartesian(xlim = c(floor(min_x), ceiling(max_x)),
                             ylim = c(floor(min_y), ceiling(max_y))) +
    scale_x_continuous(breaks = seq(from=floor(min_x), to=ceiling(max_x), by=1)) +
    scale_y_continuous(breaks = seq(from=floor(min_y), to=ceiling(max_y), by=1)) +
    ggplot2::scale_color_manual(labels=c("Activated in Males", "Inhibited in Males"),
                                values = c(RColorBrewer::brewer.pal(11, "RdYlBu")[c(1)],
                                           RColorBrewer::brewer.pal(11, "RdYlBu")[c(11)]))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("(IPA_UP) ", celltype, "_.tiff"),
                  plot = last_plot(),
                  device = "jpeg",
                  path = parent_path,
                  scale = 1,
                  width = 6,
                  height = 7,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
} 



#******************************************************************************#
# CHECK IF UPSTREAM REGULATORS OF ONE CELL TYPE IS PRODUCED BY ANOTHER CELL TYPE#
#******************************************************************************#

# cells <- c("Epithelial_Cells", "Fibroblasts", "Myeloid_Cells","B_Cells", "T_Cells" )
# 
# DEGs_df <- data.frame(matrix(ncol = 4, nrow = 0))
# colnames(DEGs_df) <- c("Gene", "cell_type", "log2FoldChange", "padj")
# upstream <- data.frame(matrix(ncol = 5, nrow = 0))
# colnames(upstream) <- c("Upstream.Regulator", "cell_type", "Expr.Log.Ratio", "Molecule.Type",
#                         "Activation.zscore")
# 
# for (cell in cells){
#   a <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "hany paper/DEGs_df_Sex_Male_vs_Female_", cell, ".xlsx"))
#   colnames(a)[1] <- "_"
#   a <- a %>% 
#     dplyr::rename(cell_type = cluster, Gene = name) %>% 
#     dplyr::filter(padj < 0.05 & log2FoldChange >= 0.58 & Gene != "None") %>%
#     dplyr::select("Gene", "cell_type", "log2FoldChange", "padj")
#   
#   DEGs_df <- rbind(DEGs_df, a)
#   
#   b <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "hany paper/DEGs_df_IPA_Upstream_", cell, ".xlsx"),
#                            startRow = 2)
#   b <- b %>% 
#     dplyr::mutate(cell_type = cell) %>% 
#     dplyr::select("Upstream.Regulator", "cell_type", "Expr.Log.Ratio", "Molecule.Type", "Activation.z-score") %>%
#     dplyr::filter(Molecule.Type == "cytokine" | Molecule.Type == "growth factor") %>%
#     dplyr::mutate_at(c("Expr.Log.Ratio", "Activation.z-score"), as.numeric)
# 
#   # Remove - from colnames as R considers them as minus when using dplyr
#   colnames(b) <- gsub("-", "", colnames(b))
#   
#   b <- b %>% dplyr::filter(Activation.zscore > 2)
#   upstream <- rbind(upstream, b)
# }
# 
# # Filter DEGs_df
# DEGS <- DEGs_df %>% dplyr::filter(Gene %in% upstream$Upstream.Regulator)




