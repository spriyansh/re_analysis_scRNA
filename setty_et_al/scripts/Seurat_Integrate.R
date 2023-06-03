#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Seurat Integrate ###########
#######################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))

# Prefix
prefixIn <- "setty_et_al/data/output/Azimuth/"
dir.create("setty_et_al/data/output/Integrated", showWarnings = F)
prefixOut <- "setty_et_al/data/output/Integrated/"

# all_file_names
all_file_names <- list.files(paste0(prefixIn))
all_file_names <- grep(all_file_names, pattern = "h5seurat", value = T)

# Make a list of seurat Object
rep_list <- list()

# Create a list
for (i in all_file_names){
    
    if (i == "rep1_Annotated_sob.h5seurat") {
        rep="rep1"
        age <- "35"
        sex <- "M"
    } else if (i == "rep2_Annotated_sob.h5seurat") {
        rep="rep2"
        age <- "28"
        sex <- "F"
    } else if (i == "rep3_Annotated_sob.h5seurat") {
        rep="rep3"
        age <- "19"
        sex <- "F"
    }
    
    #  the Processed Seurat Object
    sob <- LoadH5Seurat(file = paste0(prefixIn, i), verbose = F)
    
    sob@meta.data$sex <- sex
    sob@meta.data$age <- age
    #names(sob) <- make.unique(rep(paste(rep,age,sex, sep = "_"), length(colnames(sob))), sep = "_")
    
    # Add To list 
    rep_list <- append(rep_list, list(sob))
}

names(rep_list) <- c(paste("rep1", "35", "M", sep = "_"), 
                     paste("rep2", "28", "F", sep = "_"),
                     paste("rep3", "19", "F", sep = "_"))

# Find Vraible Features
rep_list <- lapply(rep_list, function(x){
    
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 8000)
})

# Selection Integration Feature
features <- SelectIntegrationFeatures(object.list = rep_list, nfeatures = 4000)

# Define Anchors
anchors <- FindIntegrationAnchors(object.list = rep_list, 
                                  anchor.features = features)

# Integrate the datasets
setty_integrate <- IntegrateData(anchorset = anchors)

# original unmodified data still resides in the 'RNA' assay
DefaultAssay(setty_integrate) <- "integrated"

# Run the standard workflow for visualization and clustering
setty_integrate <- ScaleData(setty_integrate, verbose = FALSE)
setty_integrate <- RunPCA(setty_integrate, npcs = 30, verbose = FALSE)
setty_integrate <- RunUMAP(setty_integrate, reduction = "pca", dims = 1:30)
setty_integrate <- FindNeighbors(setty_integrate, reduction = "pca", dims = 1:30)
setty_integrate <- FindClusters(setty_integrate, resolution = 0.5)

# # Finding Conserved Marker
# markers <- FindConservedMarkers(setty_integrate,
#                                 ident.1 = 0,
#                                 ident.2 = 5,
#                                 grouping.var = "cell.type",
#                                 verbose = FALSE)
# 
# HSC_diff <- FindMarkers(setty_integrate, ident.1 = "HSC", verbose = FALSE,
#             group.by = "cell.type")
# gn <- rownames(HSC_diff[HSC_diff$p_val <=0.05 & HSC_diff$avg_log2FC > 1,])
# 
# FeaturePlot(setty_integrate, features = gn[c(1:9)])
# 
# DotPlot(setty_integrate, features = head(rownames(markers)), 
#         cols = c("blue", "red", "green"), dot.scale = 8, split.by = "age") +
#     RotatedAxis()

# Write Seurat H5
file_name <- paste0(prefixOut, "Setty_Integrated_sob")
SaveH5Seurat(
    object = setty_integrate, filename = file_name,
    overwrite = T, verbose = FALSE
)
