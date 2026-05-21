# =========================================
# 0. 패키지 로드
# =========================================
library(Seurat)
library(SeuratDisk)

# 처음 한 번만 설치
# install.packages("SeuratDisk")

# =========================================
# 1. 데이터 로드
# =========================================
expr <- read.csv("HNSCC_expr.csv", row.names = 1)
expr_t <- t(expr)

labels <- read.csv("HNSCC_labels.csv", row.names = 1)

seurat_obj <- CreateSeuratObject(counts = expr_t)
seurat_obj$cell_type <- labels$cell_type

# =========================================
# 2. 전처리 및 임베딩 생성
# =========================================
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.9)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)


# =========================================
# 4. h5ad로 변환
# =========================================

library(Seurat)
library(Matrix)
# counts
counts <- GetAssayData(seurat_obj, slot = "counts")

# metadata
meta <- seurat_obj@meta.data

# embedding
umap <- Embeddings(seurat_obj, "umap")
pca  <- Embeddings(seurat_obj, "pca")

writeMM(counts, "matrix.mtx")
write.csv(meta, "metadata.csv")
write.csv(umap, "umap.csv")
write.csv(pca, "pca.csv")
write.csv(rownames(counts), "genes.csv", row.names = FALSE)
write.csv(colnames(counts), "cells.csv", row.names = FALSE)



write_h5ad_SeuratV5 <- function(object, file) {
  
  require(rhdf5)
  require(Seurat)
  
  ## Checks
  if (class(object) != "Seurat") stop("Please supply object.h5ad") 
  if (!stringr::str_detect(file, "h5ad")) stop("Please supply object.h5ad")
  
  ## Create file
  if (file.exists(file)) stop("File already exists.")
  h5createFile(file)
  
  # Add count matrix
  h5write(
    obj  = Matrix::as.matrix(object@assays$RNA@layers$data), 
    file = file, 
    name = "X"
  )
  
  # Add coldata (obs)
  h5write(
    obj  = cbind(cell_name = colnames(object), object@meta.data), 
    file = file, 
    name = "obs")
  
  # Add meta features (var)
  h5write(
    obj  = cbind(gene_name = rownames(object), object@assays$RNA@meta.data),
    file = file, 
    name = "var")
  
  # Add dimensional reductions (obsm)
  h5createGroup(file, "obsm")
  
  for (reduction in names(object@reductions)) {
    h5write(
      obj  = t(object@reductions[[reduction]]@cell.embeddings), 
      file = file, 
      name = paste0("obsm/X_", reduction)
    )
  }
}

write_h5ad_SeuratV5(seurat_obj2, "HNSCC_scGFT.h5ad")
























