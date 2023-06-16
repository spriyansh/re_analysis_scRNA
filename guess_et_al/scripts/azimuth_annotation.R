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
prefixIn <- "guess_et_al/data/output/p17/SeuratObjects/"
prefixIn2 <- "guess_et_al/data/output/p17/Integrated/"
prefixOut <- "guess_et_al/data/output/"

# mds genes
mds_genes <- c(
  "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
  "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
  "NF1", "CALR", "MPL", "GATA2"
)

# all_file_names
all_file_names_1 <- list.files(paste0(prefixIn))
all_file_names_2 <- list.files(paste0(prefixIn2))
all_file_names <- c(all_file_names_1, all_file_names_2)
all_file_names <- grep(all_file_names,
                       pattern = "Integrated|Processed_sob.h5seurat",
                       value = T
)

# Distribute filenames
mds17_names <- grep(all_file_names, pattern = "MDS17", value = T)
aml17_names <- grep(all_file_names, pattern = "AML17", value = T)
integrated_names <- grep(all_file_names, pattern = "Integrated", value = T)

# Make a list of filenames
rep_list <- list(
    mds17 = mds17_names,
    aml17 = aml17_names,
    integrated = integrated_names
)

# Run for every dataset
for (i in names(rep_list)) {
    
    #i = "mds17"
    
    # Get the names rep#
    rep_i_name <- rep_list[[i]]
    
    if (i == "mds17") {
        patientID <- "17"
        condition <- "MDS"
        age <- 84
        sex <- "Male"
        sob <- LoadH5Seurat(file = paste0(prefixIn, rep_i_name), verbose = F)
        
    } else if (i == "aml17") {
        patientID <- "17"
        condition <- "AML"
        age <- 84
        sex <- "Male"
        sob <- LoadH5Seurat(file = paste0(prefixIn, rep_i_name), verbose = F)
        
    }else if (i == "integrated"){
        patientID <- "17"
        condition <- "Integrated"
        age <- 84
        sex <- "Male"
        sob <- LoadH5Seurat(file = paste0(prefixIn2, rep_i_name), verbose = F)
    }
    
    # Azimuth Annotations
    bm <- RunAzimuth(
        query = sob,
        reference = paste0(
            str_sub(prefixIn, end = -26), "input/",
            "Azimuth_Human_BoneMarrow"
        )
    )
    
    # Chnage the column name
    bm@meta.data$Azimuth.Cells <- bm@meta.data$predicted.celltype.l2
    
    # Plot the Annotations
    p <- DimPlot(bm, reduction = "umap", pt.size = 1, group.by = "Azimuth.Cells") +
        ggtitle(paste0("PatientID: ", patientID, " Condition: ", condition)) +
        scale_color_hue(l = 50) + theme(legend.position = "bottom")
    outDir <- paste0(prefixOut, "p17/Azimuth/")
    dir.create(outDir, showWarnings = F)
    ggsave(p,
           filename = paste0(outDir, i, "_All_Anno_Azimuth.png"),
           dpi = 1400, limitsize = FALSE,
           width = 8, height = 8
    )
    
    # Extract Normalized Counts and New Metadata to create a CDS Object
    norm.counts <- bm@assays$RNA@data
    metadata <- as.data.frame(bm@meta.data)
    
    # Create a fresh seurat Object
    sob.fresh <- CreateSeuratObject(
        counts = bm@assays$RNA@counts,
        meta.data = metadata
    )
    # Add Normalized Counts
    sob.fresh@assays$RNA@data <- norm.counts
    
    # Add Scaled Data
    sob.fresh@assays$RNA@scale.data <- bm@assays$RNA@scale.data
    
    # Write Seurat H5
    outDir <- paste0(prefixOut, "p17/Azimuth/")
    dir.create(outDir, showWarnings = F)
    file_name <- paste0(outDir, i, "_Annotated_sob")
    SaveH5Seurat(
        object = sob.fresh, filename = file_name,
        overwrite = T, verbose = FALSE
    )
    
    # Convert to Ann Data
    #Convert(paste0(outDir, i, "_Annotated_sob.h5seurat"), dest = "h5ad", overwrite = T)
    #print(paste("Completed", i))
}
