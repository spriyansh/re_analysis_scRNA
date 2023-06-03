#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Azimuth Annotation #########
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
prefixIn <- "setty_et_al/data/input/processed/"
prefixOut <- "setty_et_al/data/output/"

# mds genes
mds_genes <- c(
  "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
  "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
  "NF1", "CALR", "MPL", "GATA2"
)

# all_file_names
all_file_names <- list.files(paste0(prefixIn))

# Distribute filenames
rep1_names <- grep(all_file_names, pattern = "rep1", value = T) # Male-35
rep2_names <- grep(all_file_names, pattern = "rep2", value = T) # Female-28
rep3_names <- grep(all_file_names, pattern = "rep3", value = T) # Female-19

# Make a list of filenames
rep_list <- list(rep1 = rep1_names, rep2 = rep2_names, rep3 = rep3_names)

# Run for every dataset
for (i in names(rep_list)) {
  if (i == "rep1") {
    individual <- "1"
    age <- "35"
    sex <- "Male"
  } else if (i == "rep2") {
    individual <- "2"
    age <- "28"
    sex <- "Female"
  } else if (i == "rep3") {
    individual <- "3"
    age <- "19"
    sex <- "Female"
  }
    
  # Get the names rep#
  rep_i_name <- rep_list[[i]]
  
  #  the Processed Seurat Object
  sob <- LoadH5Seurat(file = paste0(prefixIn, rep_i_name), verbose = F)
  
  # Azimuth Annotations
  bm <- RunAzimuth(query = sob, 
                   reference = paste0(str_sub(prefixIn, end = -11), 
                                      "Azimuth_Human_BoneMarrow")
                   )
  
  # Plot the Annotations
  p <- DimPlot(bm, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
      ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
      scale_color_hue(l = 50) + theme(legend.position = "bottom")
  
  ggsave(p,
         filename = paste0(prefixOut, i, "_All_Anno_Azimuth.png"),
         dpi = 1400, limitsize = FALSE, width = 8, height = 8
  )
  
  # Subset the Scores
  bm.score <- subset(bm, subset = predicted.celltype.l2.score >= 0.6)
  
  # Plot the Annotations
  p <- DimPlot(bm.score, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
      ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
      scale_color_hue(l = 50) + theme(legend.position = "bottom")
  
  # Save 
  ggsave(p,
         filename = paste0(prefixOut, i, "_Score_Sub_Anno_Azimuth.png"),
         dpi = 1400, limitsize = FALSE, width = 8, height = 8
  )
  
  # Subset the Cell types
  bm.cell <- subset(bm.score,
                    subset = predicted.celltype.l2 %in% c("Early Eryth",
                                                          "GMP", "HSC",
                                                          "LMPP", "Late Eryth", 
                                                          "CLP", "EMP")
                    )
  
  # Plot the Annotations
  p <- DimPlot(bm.cell, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
      ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
      scale_color_hue(l = 50) + theme(legend.position = "bottom")
  
  ggsave(p,
         filename = paste0(prefixOut, i, "_Cell_Sub_Anno_Azimuth.png"),
         dpi = 1400, limitsize = FALSE, width = 8, height = 8
  )
  
  # Plot Clusters
  p <- DimPlot(bm.cell, reduction = "umap", pt.size = 1, group.by = "seurat_clusters") +
      ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
      scale_color_hue(l = 50) + theme(legend.position = "bottom")
  
  ggsave(p,
         filename = paste0(prefixOut, i, "_Seurat_Clusters.png"),
         dpi = 1400, limitsize = FALSE, width = 8, height = 8
  )
  
  # Extract Normalized Counts and New Metadata to create a CDS Object
  norm.counts <- bm.cell@assays$RNA@data
  metadata <- as.data.frame(bm.cell@meta.data)
  metadata$cell.type <- metadata$predicted.celltype.l2
  
  # Create a fresh seurat Object
  sob.fresh <- CreateSeuratObject(counts = bm.cell@assays$RNA@counts, 
                     meta.data = metadata)
  # Add Normalized Counts
  sob.fresh@assays$RNA@data <- norm.counts
  
  # Add Scaled Data
  sob.fresh@assays$RNA@scale.data <- bm.cell@assays$RNA@scale.data
  
  # Write Seurat H5
  file_name <- paste0(prefixOut, i, "_Annotated_sob")
  SaveH5Seurat(
      object = sob.fresh, filename = file_name,
      overwrite = T, verbose = FALSE
  )
  
  # Convert to Ann Data
  Convert(paste0(prefixOut, i, "_Annotated_sob.h5seurat"), dest = "h5ad")
  
}

