#setwd("/Users/Vartika_Bisht/Individual_Project/MicrobiomeAnalystR-master/R")
#files.sources = list.files()
#sapply(files.sources, source)
#setwd("/Users/Vartika_Bisht/Individual_Project")

pacman::p_load(phyloseq, metacoder, pryr, biomformat, RColorBrewer, ggplot2, gplots, Cairo, igraph, 
               BiocParallel, randomForest, metagenomeSeq, MASS, DESeq2, vegan, RJSONIO, ggfortify, pheatmap, xtable, genefilter,
               data.table, reshape, stringr, ape, grid, gridExtra, splitstackshape, edgeR, globaltest, R.utils, viridis, ggrepel,
               ppcor)
devtools::install_github("xia-lab/MicrobiomeAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))
