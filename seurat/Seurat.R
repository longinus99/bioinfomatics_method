library(Seurat)
library(FNN)
library(entropy)
library(dplyr)


expr <- read.csv("HNSCC_expr.csv", row.names = 1)
expr_t <- t(expr)
labels <- read.csv("HNSCC_labels.csv", row.names = 1)
true_labels <- labels$cell_type

seurat_obj <- CreateSeuratObject(counts = expr_t)
head(colnames(seurat_obj))

seurat_obj$cell_type <- labels
table(seurat_obj$cell_type, useNA = "ifany")

seurat_obj1 <- NormalizeData(seurat_obj)
seurat_obj1 <- FindVariableFeatures(seurat_obj1)
seurat_obj1 <- ScaleData(seurat_obj1)
seurat_obj1 <- RunPCA(seurat_obj1)

seurat_obj1 <- RunUMAP(seurat_obj1, dims = 1:10)

DimPlot(seurat_obj1,
        reduction = "umap",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.6)
