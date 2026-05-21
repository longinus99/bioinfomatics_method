read_h5ad <- function(
    adata_path,
    embedding_key = NULL,
    gene_key = "var/gene_name",
    cell_key = "obs/_index",
    assay_name = "RNA",
    reduction_name = "scFANCL"
) {
  
  library(Seurat)
  library(rhdf5)
  library(Matrix)
  
  # expression matrix
  mat <- h5read(adata_path, "X")
  
  # names
  cells <- as.vector(
    h5read(adata_path, cell_key)
  )
  
  genes <- make.unique(
    as.vector(
      h5read(adata_path, gene_key)
    )
  )
  
  # sparse dgCMatrix
  mat <- as(
    Matrix(mat, sparse = TRUE),
    "dgCMatrix"
  )
  
  # dimnames
  rownames(mat) <- genes
  colnames(mat) <- cells
  
  # create seurat object
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    assay = assay_name
  )
  
  # embedding optional
  if (!is.null(embedding_key)) {
    
    emb <- h5read(
      adata_path,
      paste0("obsm/", embedding_key)
    )
    
    emb <- as.matrix(emb)
    
    if (nrow(emb) != length(cells)) {
      emb <- t(emb)
    }
    
    rownames(emb) <- colnames(seurat_obj)
    
    colnames(emb) <- paste0(
      reduction_name,
      "_",
      seq_len(ncol(emb))
    )
    
    seurat_obj[[reduction_name]] <- CreateDimReducObject(
      embeddings = emb,
      key = paste0(reduction_name, "_"),
      assay = assay_name
    )
  }
  
  return(seurat_obj)
}
