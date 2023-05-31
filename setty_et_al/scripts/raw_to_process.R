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
prefixIn <- "../data/input/"
prefixOut <- "../data/output/"

# mds genes
mds_genes <- c(
    "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
    "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
    "NF1", "CALR", "MPL", "GATA2"
)

# Read BioMart info
# "GO:0000087", "GO:0006260", "GO:0048285", "GO:0000279", "GO:0007059"
biomart.anno <- readRDS(paste0(prefixIn, "cell_cycle_data.mart"))
reg.out <- c(unique(biomart.anno$hgnc_symbol))

# all_file_names
all_file_names <- list.files(paste0(prefixIn))

# Distribute filenames
rep1_names <- grep(all_file_names, pattern = "rep1", value = T)
rep2_names <- grep(all_file_names, pattern = "rep2", value = T)
rep3_names <- grep(all_file_names, pattern = "rep3", value = T)

# Make a list of filenames
rep_list <- list(rep1 = rep1_names, rep2 = rep2_names, rep3 = rep3_names)

# Run for every dataset
for (i in names(rep_list)) {
    
    # Measures
    measure <- list()
    
    # Get the names rep#
    rep_i_name <- rep_list[[i]]
    rep_i_raw_name <- rep_i_name[2]
    print(rep_i_raw_name)
    rep_i_filt_name <- rep_i_name[1]
    print(rep_i_filt_name)
    
    # Load the RDS files
    raw_mat <- Read10X_h5(paste0(prefixIn, rep_i_raw_name))
    filt_mat <- Read10X_h5(paste0(prefixIn, rep_i_filt_name))
    
    # Add Measures to the list
    measure <- append(measure, list(
        raw_feature = nrow(raw_mat),
        filt_feature = nrow(filt_mat),
        raw_bc = ncol(raw_mat),
        filt_bc = ncol(filt_mat)
    ))
    # Remove Raw Mat
    raw_mat <- NULL
    
    # Renaming the features
    filt_mat@Dimnames[[1]] <- str_replace(filt_mat@Dimnames[[1]], "_", "-")
    
    # Create Seurat Object
    sob.raw <- CreateSeuratObject(
        counts = filt_mat, min.cells = 200,
        min.features = 1000, project = i
    )
    
    # Check
    match_genes <- length(rownames(sob.raw)[(rownames(sob.raw) %in% mds_genes)])
    if (match_genes == length(mds_genes)) {
        cat(paste0("\nCheck-1 (", i, "): All Genes exist\n"))
    } else if (match_genes < length(mds_genes)) {
        rem <- length(mds_genes) - match_genes
        cat(paste0("\nCheck-1 (", i, "): ", rem, " Dropped\n"))
    }
    
    # Save the raw object
    file_name <- paste0(prefixOut, i, "_raw_sob.RDS")
    saveRDS(sob.raw, file = file_name)
    
    # Calculate Percentage of Mitochondrial Read
    sob.raw[["percent.mt"]] <- PercentageFeatureSet(sob.raw, pattern = "^MT-")
    
    # Plots
    p <- VlnPlot(sob.raw, features = c("nFeature_RNA"),
                 pt.size = 1) +
        theme(title = element_text(size = 25),
            axis.text = element_text(size = 25),
            legend.position = "none"
        )
    ggsave(p, filename = paste0(prefixOut, i, "_Raw_VlnPlot_nFeature_RNA.png"),
           dpi = 1200, limitsize = FALSE)
    
    p <- VlnPlot(sob.raw, features = c("nCount_RNA"),
                 pt.size = 1) +
        theme(
            title = element_text(size = 25),
            axis.text = element_text(size = 25),
            legend.position = "none"
        )
    ggsave(p, filename = paste0(prefixOut, i, "_Raw_VlnPlot_nCount_RNA.png"),
           dpi = 1200, limitsize = FALSE)
    
    p <- VlnPlot(sob.raw, features = c("percent.mt"),
                 pt.size = 1) +
        theme(
            title = element_text(size = 25),
            axis.text = element_text(size = 25),
            legend.position = "none"
        )
    ggsave(p, filename = paste0(prefixOut, i, "_Raw_VlnPlot_percent.mt.png"),
           dpi = 1200, limitsize = FALSE)
    
    p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "percent.mt") +
        theme(
            title = element_text(size = 25), axis.text = element_text(size = 25),
            legend.position = "none"
        ) + geom_smooth(method = "lm", formula = "y~x")
    ggsave(p, filename = paste0(prefixOut, i, "_Raw_FeatureScatter_percent.mt.png"),
           dpi = 1200, limitsize = FALSE)
    
    p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        theme(
            title = element_text(size = 25), axis.text = element_text(size = 20),
            legend.position = "none"
        ) + geom_smooth(method = "lm", formula = "y~x")
    ggsave(p, filename = paste0(prefixOut, i, "_Raw_FeatureScatter_nFeature_RNA.png"),
           dpi = 1200, limitsize = FALSE)
    
    # Subset # Check Required if-else
    sob.sub <- subset(sob.raw, subset = nFeature_RNA > 200 & nCount_RNA < 30000 & percent.mt < 5)
    
    # Plots
    p <- VlnPlot(sob.sub, features = c("nFeature_RNA")) +
        theme(
            title = element_text(size = 25),
            axis.text = element_text(size = 25),
            legend.position = "none"
        )
    ggsave(p, filename = paste0(prefixOut, i, "_Sub_VlnPlot_nFeature_RNA.png"),
           dpi = 1200, limitsize = FALSE)
    
    p <- VlnPlot(sob.sub, features = c("nCount_RNA")) +
        theme(
            title = element_text(size = 25),
            axis.text = element_text(size = 25),
            legend.position = "none"
        )
    ggsave(p, filename = paste0(prefixOut, i, "_Sub_VlnPlot_nCount_RNA.png"),
           dpi = 1200, limitsize = FALSE)
    
    p <- VlnPlot(sob.sub, features = c("percent.mt")) +
        theme(
            title = element_text(size = 25),
            axis.text = element_text(size = 25),
            legend.position = "none"
        )
    ggsave(p, filename = paste0(prefixOut, i, "_Sub_VlnPlot_percent.mt.png"),
           dpi = 1200, limitsize = FALSE)
    
    p <- FeatureScatter(sob.sub, feature1 = "nCount_RNA", feature2 = "percent.mt") +
        theme(
            title = element_text(size = 25), axis.text = element_text(size = 20),
            legend.position = "none"
        ) + geom_smooth(method = "lm", formula = "y~x")
    ggsave(p, filename = paste0(prefixOut, i, "_Sub_FeatureScatter_percent.mt.png"),
           dpi = 1200, limitsize = FALSE)
    
    p <- FeatureScatter(sob.sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        theme(
            title = element_text(size = 25), axis.text = element_text(size = 20),
            legend.position = "none"
        ) + geom_smooth(method = "lm", formula = "y~x")
    ggsave(p, filename = paste0(prefixOut, i, "_Sub_FeatureScatter_nFeature_RNA.png"),
           dpi = 1200, limitsize = FALSE)
    
    
    # Check
    match_genes <- length(rownames(sob.sub)[(rownames(sob.sub) %in% mds_genes)])
    if (match_genes == length(mds_genes)) {
        cat(paste0("\nCheck-2 (", i, "): All Genes exist\n"))
    } else if (match_genes < length(mds_genes)) {
        rem <- length(mds_genes) - match_genes
        cat(paste0("\nCheck-2 (", i, "): ", rem, " Dropped\n"))
    }
    
    # Save
    file_name <- paste0(prefixOut, i, "_sub_sob.RDS")
    saveRDS(sob.sub, file = file_name)
    
    # Pre-processing
    sob.prs <- NormalizeData(sob.sub)
    sob.prs <- FindVariableFeatures(sob.prs, selection.method = "vst")
    sob.prs <- ScaleData(sob.prs, features = rownames(sob.prs))
    sob.prs <- RunPCA(sob.prs, features = VariableFeatures(sob.prs))
    
    cat(paste0("\nCheck-2 Data-pre-processed\n"))
    
    # Basic PCA
    p <- DimPlot(sob.prs, reduction = "pca")
    ggsave(p, filename = paste0(prefixOut, i, "_basic_PCA.png"),
           dpi = 1200, limitsize = FALSE)
    
    # Save
    file_name <- paste0(prefixOut, i, "_PCA_sob.RDS")
    saveRDS(sob.prs, file = file_name)
    
    # Check how many genes are present in the dataset
    indata <- rownames(sob.prs)[rownames(sob.prs) %in% reg.out]
    
    # Keep only those which are present in data
    biomart.anno <- biomart.anno[biomart.anno$hgnc_symbol %in% indata, ]
    
    # For seurat we need to divide genes into vectors
    s.genes <- unique(biomart.anno[biomart.anno$go_id %in% c("GO:0006260"), "hgnc_symbol"])
    g2m.genes <- unique(biomart.anno[biomart.anno$go_id %in% c("GO:0000087", "GO:0000279", "GO:0007059", "GO:0048285"), "hgnc_symbol"])
    
    # Cell Cycle Scores
    sob.cc <- CellCycleScoring(sob.prs, s.features = s.genes, 
                               g2m.features = g2m.genes,
                               set.ident = TRUE)
    
    # Re-Run PCA
    sob.cc <- RunPCA(sob.cc, features = c(s.genes, g2m.genes))
    
    # Save The plot
    p <- DimPlot(sob.cc)
    ggsave(p, filename = paste0(prefixOut, i, "_CCG_PCA.png"),
           dpi = 1200, limitsize = FALSE)
    
    # Scale to Regress out
    sob.cc <- ScaleData(sob.cc, vars.to.regress = c("S.Score", "G2M.Score"), 
                        features = rownames(sob.cc))
    
    # Re-Run PCA
    sob.cc <- RunPCA(sob.cc, features = VariableFeatures(sob.cc), nfeatures.print = 10)
    
    # Save The plot
    p <- DimPlot(sob.cc)
    ggsave(p, filename = paste0(prefixOut, i, "_CCC_PCA.png"),
           dpi = 1200, limitsize = FALSE)
    
    # Save
    file_name <- paste0(prefixOut, i, "_CCC_sob.RDS")
    saveRDS(sob.cc, file = file_name)
    
    # More Pre-processing
    sob.p <- FindNeighbors(sob.cc)
    sob.p <- FindClusters(sob.p)
    sob.p <- RunUMAP(sob.p, dims = 1:30)
    
    file_name <- paste0(prefixOut, i, "_Processed_sob.RDS")
    saveRDS(sob.p, file = file_name)
    
    # Save the UMAP with the cluster
    p <- DimPlot(sob.p, reduction = 'umap')
    ggsave(p, filename = paste0(prefixOut, i, "_C_UMAP.png"),
           dpi = 1200, limitsize = FALSE)
    
    #---Azimuth
    azi.anno <- RunAzimuth(query = file_name,
                           reference = paste0(prefixIn, "Azimuth_Human_BoneMarrow"))
    
    # Save
    file_name <- paste0(prefixOut, i, "_Azimuth_sob.RDS")
    saveRDS(azi.anno, file = file_name)
    
    # Add Measures to the list
    measure <- append(measure, list(
        processed_feature = nrow(azi.anno),
        processed_bc = ncol(azi.anno)
    ))
    
    # Save
    write.table(t(as.data.frame(measure)),
                paste0(prefixOut, i,"_Measure.tsv", sep = "\t"),
                col.names = NA)
    
    cat(paste0("\nCompleted for ", i, "\n"))
}
