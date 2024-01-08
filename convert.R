library(Seurat)

rhap <- readRDS(file = "Sample_Seurat.rds")
SaveH5Seurat(rhap, filename = "Sample_Seurat.h5Seurat")
Convert("Sample_Seurat.h5Seurat", dest = "h5ad")