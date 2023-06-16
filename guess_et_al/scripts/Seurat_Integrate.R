#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Seurat Integrate p17 #######
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

bar_theme <- theme_classic()+ 
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.5)),
          axis.text.y = element_text(angle = 45, hjust = 1, size = rel(1.5)),
          axis.title = element_text(size = rel(2)),
          legend.text = element_text(size = rel(3)),
          legend.title = element_text(size = rel(4)))

mds_genes <- c(
    "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
    "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
    "NF1", "CALR", "MPL", "GATA2"
)

# Prefix
prefixIn <- "guess_et_al/data/output/p17/SingleR/"
dir.create("guess_et_al/data/output/p17/Integrated/", showWarnings = F)
prefixOut <- "guess_et_al/data/output/p17/Integrated/"

# all_file_names
all_file_names <- list.files(paste0(prefixIn))
all_file_names <- grep(all_file_names, pattern = "aml|mds", value = T)
all_file_names <- grep(all_file_names, pattern = "seurat", value = T)

# Make a list of seurat Object
rep_list <- list()

# Create a list
for (i in all_file_names) {
    
    #i = "AML17_Processed_sob.h5seurat"
    
    if (i == "aml17_Annotated_sob.h5seurat") {
        patientID <- "17"
        condition <- "AML"
        age <- 84
        sex <- "Male"
    } else if (i == "mds17_Annotated_sob.h5seurat") {
        patientID <- "17"
        condition <- "MDS"
        age <- 84
        sex <- "Male"
    }
    
    #  the Processed Seurat Object
    sob <- LoadH5Seurat(file = paste0(prefixIn, i), verbose = F)
    
    # Add meta information
    sob@meta.data$condition <- condition
    sob@meta.data$sex <- sex
    sob@meta.data$age <- age
    
    # Add To list
    rep_list <- append(rep_list, list(sob))
}

names(rep_list) <- c("AML", "MDS")

# Extract MetaData from both files
meta <- rbind(rep_list[["AML"]]@meta.data, rep_list[["MDS"]]@meta.data)
meta <- meta %>% select(c(orig.ident, Azimuth.Cells))

df_count <- as.data.frame(table(meta$orig.ident, meta$Azimuth.Cells))
names(df_count) <- c("condition", "Azimuth.Cells", "count")
bar <- ggplot(df_count, aes(fill=Azimuth.Cells, y=count, x=condition)) + 
    geom_bar(position="stack", stat="identity") + theme_minimal()

outDir <- paste0(prefixOut, "Bars/")
dir.create(outDir, showWarnings = F)
ggsave(bar, filename = paste0(outDir,"CellType_Bar.png"),
       dpi = 500, limitsize = FALSE, width = 5
)

metaSub <- meta %>% filter(Azimuth.Cells %in% c("HSC", "GMP", "LMPP"))
df_count <- as.data.frame(table(metaSub$orig.ident, metaSub$Azimuth.Cells))
names(df_count) <- c("condition", "Azimuth.Cells", "count")
bar <- ggplot(df_count, aes(fill=Azimuth.Cells, y=count, x=condition)) + 
    geom_bar(position="stack", stat="identity") + theme_minimal()

outDir <- paste0(prefixOut, "Bars/")
dir.create(outDir, showWarnings = F)
ggsave(bar, filename = paste0(outDir,"CellTypeSub_Bar.png"),
       dpi = 500, limitsize = FALSE, width = 5
)

aml <- rep_list[["AML"]]
mds <- rep_list[["MDS"]]

# Subset on the expression level of a gene/feature
mds <- subset(x = mds, subset = Azimuth.Cells %in% c("HSC", "LMPP", "GMP"))
aml <- subset(x = aml, subset = Azimuth.Cells %in% c("HSC", "LMPP", "GMP"))

list_cells <- list(AML = aml,
                   MDS = mds)

# Find Variable Features
list_cells <- lapply(list_cells, function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000,
                              verbose = F)
})

# Selection Integration Feature
features <- SelectIntegrationFeatures(object.list = list_cells, nfeatures = 2000)

# Define Anchors
anchors <- FindIntegrationAnchors(
    object.list = list_cells,
    anchor.features = features
)

# Integrate the datasets
guess_integrate <- IntegrateData(anchorset = anchors)

# original unmodified data still resides in the 'RNA' assay
DefaultAssay(guess_integrate) <- "integrated"

# Run the standard workflow for visualization and clustering
guess_integrate <- ScaleData(guess_integrate, verbose = FALSE)
guess_integrate <- RunPCA(guess_integrate, npcs = 50, verbose = FALSE)
guess_integrate <- RunUMAP(guess_integrate, reduction = "pca", dims = 1:50)
guess_integrate <- FindNeighbors(guess_integrate, reduction = "pca", dims = 1:50)
guess_integrate <- FindClusters(guess_integrate, resolution = 0.5)

# Explore Integration
p <- DimPlot(guess_integrate, group.by = "Azimuth.Cells", pt.size = 0.5) +
    ggtitle(paste0(
        "PatientID: ", patientID, " Integrated Cells", " Age: ", age,
        " Sex: ", sex
    ))+ scale_color_brewer(palette = "Dark2") +
    theme(title = element_text(size = 10))
outDir <- paste0(prefixOut, "UMAP/")
dir.create(outDir, showWarnings = F)
ggsave(p, filename = paste0(outDir, str_remove(i, "_Processed_sob.h5seurat"), "_Integrated_UMAP.png"),
       dpi = 1200, limitsize = FALSE, width = 5
)

sob.fresh <- CreateSeuratObject(counts = guess_integrate@assays$integrated@data,
                                meta.data = guess_integrate@meta.data)

# Write Seurat H5
file_name <- paste0(prefixOut, "Guess_Integrated_sob")
SaveH5Seurat(
    object = sob.fresh, filename = file_name,
    overwrite = T, verbose = FALSE
)

# Show compartments
umap_coords <- ggplot_build(p)$data[[1]]
umap_coords$alpha <- sob.fresh@meta.data$orig.ident
umap_coords$cellType <- sob.fresh@meta.data$Azimuth.Cells

a <- ggplot(umap_coords,
       aes(x = x, y = y, color = cellType, alpha = alpha)) +
    scale_alpha_manual(values = c(MDS17 = 1, AML17 = 0.2))+
    geom_point(size = 1) + theme_minimal()+ scale_color_brewer(palette = "Dark2")+
    theme(legend.position = "none") + ggtitle("MDS State")

b <- ggplot(umap_coords,
            aes(x = x, y = y, color = cellType, alpha = alpha)) +
    scale_alpha_manual(values = c(MDS17 = 0.2, AML17 = 1))+
    geom_point(size = 1) + theme_minimal()+ scale_color_brewer(palette = "Dark2")+
    theme(legend.position = "none") + ggtitle("AML State")
p <- a +b


ggsave(p, filename = paste0(outDir, "P17_Integrated_UMAP.png"),
       dpi = 1200, limitsize = FALSE, width = 8
)
