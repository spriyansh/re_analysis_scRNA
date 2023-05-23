# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))

# Prefix
prefixIn <- "../input/"
prefixOut <- "../output/"

# mds genes
mds_genes <- c(
    "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
    "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
    "NF1", "CALR", "MPL", "GATA2"
)

# all_file_names
all_file_names <- list.files(paste0(prefixOut))

# Distribute filenames
azi_obj_name <- grep(all_file_names, pattern = "Azimuth", value = T)

# Run for every dataset
for (i in azi_obj_name) {
    
    # Get the Replicate
    rep_split <- str_split(i, "_")
    names(rep_split) <- rep_split[[1]][1]
    
    # Read the dataset
    azi_obj <- readRDS(paste0(prefixOut, i))
    prs_obj <- readRDS(paste0(prefixOut,rep_split[[1]][1],"_Processed_sob.RDS"))
    
    # Add Labels
    azi_meta <- as.data.frame(azi_obj@meta.data)
    prs_obj$az_labels <- azi_meta$predicted.celltype.l2
    prs_obj$az_labels_score <- azi_meta$predicted.celltype.l2.score
    
    # Plot UMAP
    png(paste0(prefixOut, rep_split[[1]][1], "_UMAP_AZ_Cells.png"), width = 500, height = 500)
    p <- DimPlot(prs_obj, reduction = "umap", 
            group.by = "az_labels", pt.size = 2)
    print(p)
    dev.off()
    
    # Distribution
    png(paste0(prefixOut, rep_split[[1]][1], "_HIST_AZ_Score.png"), width = 500, height = 500)
    hist(prs_obj$az_labels_score)
    abline(v=0.6, col = "red")
    dev.off()
    
    # Subset
    prs_obj_sub <- subset(prs_obj, subset = az_labels_score > 0.6)

    # Subset based on score
    png(paste0(prefixOut, rep_split[[1]][1], "_UMAP_AZ_SubScore.png"), width = 500, height = 500)
    p <- DimPlot(prs_obj_sub, reduction = "umap", 
            group.by = "az_labels", pt.size = 2)
    print(p)
    dev.off()
    
    # Selection Based on Cell type
    cell.types <- as.data.frame(t(table(prs_obj_sub$az_labels)))
    cell.types <- cell.types[,-1]
    colnames(cell.types) <- c("cell_type", "number")
    
    # Select the required cells
    keep.cells <- c("CLP", "Early Eryth",
                    "EMP", "GMP", "HSC",
                    "Late Eryth", "LMPP")
    
    # Subset
    prs_obj_sub <- subset(prs_obj_sub, subset = az_labels %in% keep.cells)
    
    # Save Plot Again
    png(paste0(prefixOut, rep_split[[1]][1], "_UMAP_AZ_SubCell.png"), width = 500, height = 500)
    p <- DimPlot(prs_obj_sub, reduction = "umap", 
            group.by = "az_labels", pt.size = 2)
    print(p)
    dev.off()
    
    # Save Data 
    rds_name <- paste0(prefixOut, rep_split[[1]][1], "_final_obj.RDS")
    saveRDS(prs_obj_sub, rds_name)
    
    cat(paste0("\nCompleted for ", rep_split[[1]][1], "\n"))
}
