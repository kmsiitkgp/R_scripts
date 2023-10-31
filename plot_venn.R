#!/usr/bin/env Rscript

.libPaths("/hpc/home/kailasamms/NGSTools/R_packages")
.libPaths()
chooseCRANmirror(ind=1)

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

# if (!require("BiocManager", quietly = TRUE)) 
#   install.packages("BiocManager")       # Needed to download BioManager packages
# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")           # Needed to download github packages

# # Data analysis packages
# BiocManager::install("TCGAbiolinks")        # Needed for TCGA data analysis
# BiocManager::install("ensembldb")           # Needed for annotating genes
# BiocManager::install("AnnotationHub")       # Needed for annotating genes
# BiocManager::install("affy")                # Needed for micro-array analysis
# BiocManager::install("lumi")                # Needed for micro-array analysis
# BiocManager::install("illuminaHumanv4.db")  # Needed for micro-array analysis
# BiocManager::install("limma")               # Needed for micro-array analysis
# BiocManager::install("ChIPQC")              # Needed for ChIP analysis
# BiocManager::install("clusterProfiler")     # Needed for GSEA analysis
# BiocManager::install("fgsea")               # Needed for GSEA analysis
# BiocManager::install("enrichplot")          # Needed for GSEA analysis
# BiocManager::install("pathview")            # Needed for GSEA analysis  
# install.packages("msigdbr")                 # Needed for GSEA analysis
# BiocManager::install("DESeq2")              # Needed for Differential analysis
# BiocManager::install("apeglm")              # Needed for Differential analysis
# install.packages("ashr")                    # Needed for Differential analysis

# # Data wrangling packages
# install.packages("openxlsx")            # Needed for reading, writing xlsx files
# install.packages("dplyr")               # Needed for data wrangling
# install.packages("tibble")              # Needed for advanced data wrangling
# install.packages("stringr")             # Needed for advanced data wrangling
# install.packages("purrr")               # Needed for advanced data wrangling
# install.packages("Matrix")              # Needed for numerical data wrangling

# # Graph plotting packages
# install.packages("ggplot2")             # Needed for making graphs
# install.packages("cowplot")             # Needed for merging multiple graphs
# install.packages("viridis")             # Needed for nice graph coloring
# install.packages("RColorBrewer")        # Needed for nice graph coloring
# install.packages("ggrepel")             # Needed for making graphs prettier
# install.packages("ggpubr")              # Needed for adding p values to graphs
# install.packages("extrafont")           # Needed for using better fonts

# # Specialized Graph plotting packages
# BiocManager::install("ComplexHeatmap")  # Needed for making heatmap
# install.packages("pheatmap")            # Needed for making heatmap 
# install.packages("VennDiagram")         # Needed for making Venn diagram
# install.packages("survival")            # Needed for making survival curves
# install.packages("survminer")           # Needed for making survival curves

# # Single cell analysis packages
# BiocManager::install('multtest')        # Needed for installing metap
# install.packages('metap')               # Needed for using FindConservedMarkers
# install.packages("Seurat")              # Needed for single cell analysis
# install.packages("devtools")
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
library("TCGAbiolinks")
library("ensembldb")
library("AnnotationHub")
library("affy")
library("lumi")
library("illuminaHumanv4.db")
library("limma")
library("ChIPQC")
library("clusterProfiler")
library("fgsea")
library("enrichplot")
library("pathview")
library("msigdbr")
library("DESeq2")
library("apeglm")
library("ashr")

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")
library("stringr")
library("purrr")
library("Matrix")

# Graph plotting packages
library("ggplot2")
library("cowplot")
library("viridis")
library("RColorBrewer")
library("ggrepel")
library("ggpubr")
library("extrafont")

# Specialized Graph plotting packages
library("ComplexHeatmap")
library("pheatmap")
library("VennDiagram")
library("survival")
library("survminer")

# Single cell analysis packages
library("metap")
library("Seurat")
library("DoubletFinder")

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

#REQUIREMENTS:
#(i) ONE xlsx file in parent_path.
#(ii) The data to plotted MUST be in the first sheet of excel file
#(iii) The data MUST have headers corresponding to how you want to label the 
# circles in Venn Diagram
#(iv) Maximum of 5 data sets can be plotted

parent_path <- "C:/Users/KailasammS/Desktop/venn/"

# Choose one of the below three palettes with upto 5 colorblind colors:
palette1 <- c("#ffb6db","#db6d00","#490092","#006ddb","#004949") # regular
palette1 <- c("#924900","#009292","#ff6db6","#24ff24","#b66dff")
palette1 <- c("#6db6ff","#b6dbff","#920000","#000000","#ffff6d")
palette1 <- RColorBrewer::brewer.pal(11, "RdYlBu")
#palette1 <- c("#00008C", "#E75480", "#db6d00","#490092", "#313695", "#A50026") #
#pie(rep(1,length(palette1)), col=palette1)

# Give plot title
plot_title <- ""

#******************************************************************************#
#                              PLOT VENN DIAGRAM                               #
#******************************************************************************#

# Use this website for 3 set data if this script's venn diagram isnt pretty
# http://www.ehbio.com/test/venn/#/

# Import data
file <- list.files(path = parent_path)
file <- file[grepl(pattern=".xlsx", x=file)]
data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, file))

# Number of columns being plotted
ncol <- ncol(data)

# Fix column names
colnames(data) <- base::gsub("_", " ", colnames(data))
colnames(data) <- base::gsub("\\.", " ", colnames(data))

# Declare cat.pos and cat.dis
if (ncol == 4){
  pos <- c(330, 15, 330, 15)
  dist <- c(0.27, 0.25, 0.15, 0.13)
  cex = 2
  palette1 <- c("#C8E7F5", "#00008C", "#F6D2E0", "#E75480")        # for 5B, 6B
} else if (ncol == 3){
    pos <- c(0, 0, 180)
    dist <- c(0.1, 0.1, 0.1)
    cex = 2
    palette1 <- c("#C8E7F5", "#F6D2E0", "#db6d00")    # c("#00008C", "#E75480", "#db6d00") # for 7A
} else if (ncol == 2){
  pos <- c(0, 0)
  dist <- c(0.05, 0.05)
  cex = 2.75
  palette1 <- c("#C8E7F5", "#db6d00")                              # for 7B
} else if (ncol == 1){
  pos <- c(0)
  dist <- c(0.1)
  cex = 2.75
  palette1 <- c("#F6D2E0")                                         # for 7B
}

# Create a dataframe where we store the wrapped column names
annotation <- data.frame(colnames(data))
colnames(annotation) <- c("Labels")
annotation <- annotation %>% 
  dplyr::mutate(Labels = stringr::str_wrap(Labels, 6))  # Set 10 for #5C. Normally, 15

# Convert the data frame to a named list
genes <- base::vector(mode = "list", length = ncol(data))
names(genes) <- annotation$Labels

for (i in 1:ncol(data)){
  # remove NA values and create a list of genes for each label
  genes[[i]] <- data[!is.na(data[i]),i]
}

# Plot the venn diagram
VennDiagram::venn.diagram(x = genes,
                          main = plot_title, 
                          category.names = annotation$Labels,
                          filename = paste0(parent_path, "Venn_Diagram.tiff"),
                          output = TRUE,
                          scaled =FALSE,
                          imagetype = "tiff",
                          height = 11, 
                          width = 11,
                          units = "in",
                          resolution = 600,
                          compression = "lzw",
                          margin = 0.3,    #amount of white space around Venn Diagram in grid units
                          
                          # Formatting the shapes of venn diagram
                          lwd = 1.5,                 #thickness of line
                          lty = 1,                   #type of line
                          col = "black",             #color of line
                          
                          # Formatting numbers inside venn diagram
                          cex = cex,                 #font size (2 or 2.75)
                          fontface = "bold",         #font style
                          fontfamily = "sans",       #font type
                          
                          # Formatting title of venn diagram
                          main.cex = 2,              #font size
                          main.fontface = "bold",    #font style
                          main.fontfamily = "sans",  #font type
                          main.col = "black",        #font color
                          
                          # Formatting category of venn diagram
                          cat.cex = 2,               #font size
                          cat.fontface = "bold",     #font style
                          cat.fontfamily = "sans",   #font type
                          cat.col = palette1[c(1:ncol)],  #"black",  #c("#00008C", "#00008C", "#E75480", "#E75480"),
                          
                          # Formatting colors of venn diagram
                          fill = palette1[c(1:ncol)],
                          alpha = rep(0.5, ncol), #0.5=50% transparency, 1=0% transparency
                          #cat.default.pos = "outer",    
                          
                          cat.pos = pos,    
                          cat.dist = dist, 
                          disable.logging = FALSE,
                          ext.text = TRUE)

#******************************************************************************#
#                          SAVE THE OVERLAPPING GENES                          #
#******************************************************************************#

# Save the list of overlapping genes. NOTE: You need to manually figure out
# which genes belong to which overlap based on number of genes overlapping
overlap <- VennDiagram::calculate.overlap(x = genes)

# Identify maximum number of genes present in any overlap
max = max(lengths(overlap))

# Create an dataframe of size length(overlap), max with NAs
results = data.frame(matrix("", nrow = max, ncol = length(overlap)))

rownames(results) <- paste0("Gene#", seq(max))
colnames(results) <- paste0("Intersection#", seq(length(overlap)))

# Populate the dataframe with gene names
for (i in 1:length(overlap)){
  if (length(overlap[[i]]) > 0){
    for (j in 1:length(overlap[[i]])){
      results[[j,i]] <- overlap[[i]][j]
    }
  }
}

# Save the results
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, sheetName = "venn_input")
openxlsx::writeData(wb, sheet = "venn_input", x = data, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "venn_output")
openxlsx::writeData(wb, sheet = "venn_output", x = results, rowNames = TRUE)

openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Overlap.xlsx"), overwrite = TRUE,  returnValue = FALSE)
