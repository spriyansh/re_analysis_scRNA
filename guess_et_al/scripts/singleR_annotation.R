#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: SingleR Annotation #########
#######################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))

# Prefix
prefixIn <- "guess_et_al/data/output/p17/Azimuth/"
prefixOut <- "guess_et_al/data/output/"

# mds genes
mds_genes <- c(
    "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
    "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
    "NF1", "CALR", "MPL", "GATA2"
)

# all_file_names
all_file_names <- list.files(paste0(prefixIn))
all_file_names <- grep(all_file_names,
                       pattern = ".h5seurat",
                       value = T
)

# Cell Dex Reference
cell_dex_NovershternHematopoieticData <- readRDS("guess_et_al/data/input/CellDex_NovershternHematopoieticData.RDS")

# Distribute filenames
mds17_names <- grep(all_file_names, pattern = "mds17", value = T)
aml17_names <- grep(all_file_names, pattern = "aml17", value = T)
integrated_names <- grep(all_file_names, pattern = "integrated", value = T)

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
        
    } else if (i == "aml17") {
        patientID <- "17"
        condition <- "AML"
        age <- 84
        sex <- "Male"
        
    }else if (i == "integrated"){
        patientID <- "17"
        condition <- "Integrated"
        age <- 84
        sex <- "Male"
    }
    
    # Load Seurat Object
    sob <- LoadH5Seurat(file = paste0(prefixIn, rep_i_name), verbose = F)
    
    # SingleR Labels
    SingleR_labels <- SingleR(test = GetAssayData(sob, slot = "data"),
                              ref = cell_dex_NovershternHematopoieticData,
                              labels = cell_dex_NovershternHematopoieticData$label.main)
    
    # Chnage the column name
    sob@meta.data$SingleR.Cells <- SingleR_labels$labels
    
    # Calculate low dimensions
    sob.prs <- NormalizeData(sob, verbose = F)
    sob.prs <- FindVariableFeatures(sob.prs, selection.method = "vst", verbose = F,
                                    nfeatures = 6000)
    sob.prs <- ScaleData(sob.prs, features = rownames(sob.prs), verbose = F)
    sob.prs <- RunPCA(sob.prs,
                      features = VariableFeatures(sob.prs),
                      verbose = F, ndims.print = 0, nfeatures.print = 0
    )
    sob.p <- FindNeighbors(sob.prs, verbose = F)
    sob.p <- FindClusters(sob.p, verbose = F)
    sob <- RunUMAP(sob.p, dims = 1:10, verbose = F)
    
    # Plot the Annotations
    p <- DimPlot(sob, reduction = "umap", pt.size = 1, group.by = "SingleR.Cells") +
        ggtitle(paste0("PatientID: ", patientID, " Condition: ", condition)) +
        scale_color_hue(l = 50) + theme(legend.position = "bottom")
    outDir <- paste0(prefixOut, "p17/SingleR/")
    dir.create(outDir, showWarnings = F)
    ggsave(p,
           filename = paste0(outDir, i, "_All_Anno_SingleR.png"),
           dpi = 1400, limitsize = FALSE,
           width = 8, height = 8
    )
    
    # Write Seurat H5
    outDir <- paste0(prefixOut, "p17/SingleR/")
    dir.create(outDir, showWarnings = F)
    file_name <- paste0(outDir, i, "_Annotated_sob")
    SaveH5Seurat(
        object = sob, filename = file_name,
        overwrite = T, verbose = FALSE
    )
}
