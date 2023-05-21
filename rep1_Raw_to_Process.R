# Load the Required Libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))








## Rename File












# Load the raw data
s.filt.mtx <- Read10X_h5("raw_data/filtered_feature_bc_matrix.h5")

# Renaming the features
s.filt.mtx@Dimnames[[1]] <- stringr::str_replace(s.filt.mtx@Dimnames[[1]], "_", "-")

# Create Seurat Object
rep1.raw <- CreateSeuratObject(counts = s.filt.mtx, min.cells = 200,
                               min.features = 1000, project = "rep1") 
saveRDS(rep1.sob, "processed_data/rep1.raw.RDS")

# Calculate Percentage of Mitochondrial Read
rep1.raw[["percent.mt"]] <- PercentageFeatureSet(rep1.raw, pattern = "^MT-")

# Subset
rep1.sub <- subset(rep1.raw, subset = nFeature_RNA > 200 & nCount_RNA < 45000 & percent.mt < 10)

# Normalize
rep1.sub <- NormalizeData(rep1.sub)

# Find HVGs
rep1.sub <- FindVariableFeatures(rep1.sub, selection.method = "vst")

# Scale Data
rep1.sub <- ScaleData(rep1.sub, features = rownames(rep1.sub))

# Run PCA
rep1.sub.pca <- RunPCA(rep1.sub, features = VariableFeatures(rep1.sub))
saveRDS(rep1.sub.pca, "processed_data/rep1.sub.pca.RDS")

# Read Mart Data
biomart.anno <- readRDS("raw_data/cell_cycle_data.mart")

# Genes that are present in the dataset
indata <- rownames(rep1.sub.pca)[rownames(rep1.sub.pca) %in% unique(biomart.anno$hgnc_symbol)]
biomart.anno <- biomart.anno[biomart.anno$hgnc_symbol %in% indata,]

# For seurat we need to divide genes into vectors
s.genes <- unique(biomart.anno[biomart.anno$go_id %in% c("GO:0006260"), "hgnc_symbol"])
g2m.genes <- unique(biomart.anno[biomart.anno$go_id %in% c("GO:0000087", "GO:0000279", "GO:0007059", "GO:0048285"), "hgnc_symbol"])

# Cell Cycle Scores
rep1.sub.pca.cc <- CellCycleScoring(rep1.sub.pca, 
                                     s.features = s.genes,
                                     g2m.features = g2m.genes,
                                     set.ident = TRUE)

# Run PCA again
rep1.sub.pca.cc <- RunPCA(rep1.sub.pca.cc,
                           features = c(s.genes, g2m.genes))

saveRDS(rep1.sub.pca.cc, "processed_data/rep1.sub.pca.cc.RDS")

# Scale Data
rep1.sub.pca.cc.removed <- ScaleData(rep1.sub.pca.cc, 
                              vars.to.regress = c("S.Score", "G2M.Score"), 
                              features = rownames(rep1.sub.pca.cc))

# Run PCA again to check for 
rep1.sub.pca.cc.removed <- RunPCA(rep1.sub.pca.cc.removed, features = VariableFeatures(rep1.sub.pca.cc.removed), 
                                  nfeatures.print = 10)

# SaveRDS
saveRDS(rep1.sub.pca.cc.removed. "processed_data/rep1.sub.pca.cc.removed.RDS")