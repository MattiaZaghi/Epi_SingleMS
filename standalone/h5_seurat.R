# Transform Anndata to Seurat
.libPaths("/home/mattia/miniconda3/envs/scanpy/lib/R/library/")
library("reticulate")
library("Seurat")
use_condaenv("scanpy")

sc <- import("scanpy")

# Load data from h5ad
adata <- sc$read("/date/gcb/gcb_MZ/multiNanoCT/dev_all.loom")
exprs <- t(adata$layers['spliced'] + adata$layers['unspliced'])
colnames(exprs) <- adata$obs_names$to_list()
rownames(exprs) <- adata$var_names$to_list()

# Create the Seurat object
seurat <- CreateSeuratObject(counts=exprs)

# Set normalized data
normalized <- t(adata$X)
colnames(normalized) <- adata$obs_names$to_list()
rownames(normalized) <- adata$var_names$to_list()

# Set the expression assay
seurat <- SetAssayData(seurat, "data", normalized)

# Add observation metadata
seurat <- AddMetaData(seurat, adata$obs)

# Add embedding
embedding <- adata$obsm["X_umap"]
rownames(embedding) <- adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")

# Add embedding
embedding <- adata$obsm["X_pca"]
rownames(embedding) <- adata$obs_names$to_list()
seurat[["pca"]] <- CreateDimReducObject(embedding, key = "pca_")

# Add embedding
embedding <- adata$obsm["X_pca_harmony"]
rownames(embedding) <- adata$obs_names$to_list()
seurat[["pca_harmony"]] <- CreateDimReducObject(embedding, key = "harmony")

# Export
saveRDS(file = "hca_dev_opcs_raw.rds", seurat)
test <- readRDS(file = "hca_dev_opcs_raw.rds")
