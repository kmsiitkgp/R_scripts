#!/usr/bin/env Rscript

.libPaths(new = "/hpc/home/kailasamms/R/x86_64-pc-linux-gnu-library/4.1")
.libPaths()
chooseCRANmirror(ind=1)

# NOTE: Install all packages from submit node NOT computing node. 
# NOTE: These packages CANNOT be installed from within R. Use conda. They are
# needed for ktplots()
conda install -c conda-forge r-textshaping
conda install -c conda-forge r-ragg
conda install -c conda-forge r-gert
conda install -c conda-forge r-rjags

# NOTE: Some packages like Seurat-disk, ScopeLoomR need hdf5r package which in 
# turn needs HDF5 library files to install as well as even work. Similarly, 
# metap needs fftw3. So, 
# (i) qrsh -> conda activate R -> module load hdf5/1.8.18 -> R -> install... 
# (ii) qrsh -> conda activate R -> conda install -c conda-forge fftw -> R -> install...

# NOTE: If you get the following error, 
# "configure: error: "libxml not found"
# "ERROR: configuration failed for package ‘XML’"
# Exit R, then activate conda environment and type "conda install r-xml"
# IMPORTANT: DO NOT UPDATE XML even if outdated. 

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

#"illuminaHumanv4.db", "pathview",

# Non-CRAN repository managing packages
install.packages(pkgs = c("BiocManager", "remotes"),
                 repos ='http://cran.us.r-project.org',
                 INSTALL_opts = '--no-lock')

# Data analysis packages
# BiocNeighbors is key dependecy for CellChat package
BiocManager::install(pkgs = c("TCGAbiolinks", "ensembldb", "AnnotationHub",
                              "affy", "lumi", "ChIPQC", "fgsea", "enrichplot",
                              "clusterProfiler", "DESeq2", "progeny",  
                              "dorothea", "viper", "org.Hs.eg.db", 
                              "org.Mm.eg.db", "BiocNeighbors", "infercnv"),
                     INSTALL_opts = '--no-lock')

install.packages(pkgs = c("msigdbr", "wordcloud"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Data wrangling packages
install.packages(pkgs = c("openxlsx", "dplyr", "tibble", "stringr", "purrr"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Graph plotting packages
install.packages(pkgs = c("ggplot2", "cowplot", "viridis", "RColorBrewer", "colorspace"),
                 repos = 'http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Specialized Graph plotting packages
install.packages(pkgs = c("pheatmap", "ggridges", "VennDiagram", "survival", "survminer", "ggbeeswarm"),
                 repos = 'http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Interface with Python packages
install.packages(pkgs = c("reticulate"),
                 repos = 'http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Single cell analysis packages
# arrow package is key dependency for SCENIC
# Install using compute node. It might fail using submit node.
Sys.setenv("ARROW_R_DEV"=TRUE, "LIBARROW_BINARY"=FALSE,
           "ARROW_WITH_ZSTD"="ON", "ARROW_DEPENDENCY_SOURCE"="BUNDLED")
install.packages(pkgs = c("arrow"),
                 repos = 'http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

install.packages(pkgs = c("Seurat","metap"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')
remotes::install_github(repo = c("hhoeflin/hdf5r",
                                 "mojaveazure/seurat-disk",
                                 "chris-mcginnis-ucsf/DoubletFinder",
                                 "neurorestore/Augur",
                                 "aertslab/SCENIC",
                                 "aertslab/SCopeLoomR",
                                 "zktuong/ktplots",
                                 "sqjin/CellChat"),
                        INSTALL_opts = '--no-lock')

# old.packages(lib.loc = "/hpc/home/kailasamms/R/x86_64-pc-linux-gnu-library/4.1")
# update.packages(#lib.loc = "/hpc/home/kailasamms/R/x86_64-pc-linux-gnu-library/4.1",
#                 INSTALL_opts = '--no-lock',
#                 ask = FALSE)

pkgs <- c("BiocManager", "remotes", "TCGAbiolinks", "ensembldb", "AnnotationHub",
          "affy", "lumi", "ChIPQC", "fgsea", "enrichplot", "clusterProfiler",
          "org.Hs.eg.db", "org.Mm.eg.db", "DESeq2", "progeny", "dorothea", 
          "viper", "msigdbr", "ashr", "openxlsx", 
          "dplyr", "tibble", "stringr", "purrr", "ggplot2", "cowplot", 
          "viridis", "RColorBrewer", "pheatmap", "ggridges", "VennDiagram", 
          "survival", "survminer", "ggbeeswarm", "reticulate", "arrow", 
          "multtest", "RcisTarget", "AUCell", "Seurat", "metap", "hdf5r", 
          "SeuratDisk", "DoubletFinder", "Augur", "SCENIC", "SCopeLoomR", 
          "ktplots", "CellChat", "wordcloud", "infercnv")

# Display packages that couldn't be installed
cat("\nNOT Installed:", pkgs[!(pkgs %in% installed.packages()[,1])], "\n")