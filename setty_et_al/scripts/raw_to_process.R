#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Pre-Process and CC Remove ##
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
prefixIn <- "setty_et_al/data/input/"
prefixOut <- "setty_et_al/data/output/"

# mds genes
mds_genes <- c(
  "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
  "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
  "NF1", "CALR", "MPL", "GATA2"
)

# Read BioMart info
biomart.anno <- readRDS(paste0(prefixIn, "cell_cycle_data.mart"))
reg.out <- c(unique(biomart.anno$SYMBOL))

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
  # i <- "rep1"

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

  # Measures
  measure <- list()

  # Get the names rep#
  rep_i_name <- rep_list[[i]]
  rep_i_raw_name <- rep_i_name[2]
  rep_i_filt_name <- rep_i_name[1]

  # Load the RDS files
  raw_mat <- Read10X_h5(paste0(prefixIn, rep_i_raw_name))
  filt_mat <- Read10X_h5(paste0(prefixIn, rep_i_filt_name))

  cat(paste0("\nCheck-1 (", i, "): ", "Loaded Files\n"))

  # Add Measures to the list
  measure <- append(measure, list(
    raw_feature = nrow(raw_mat),
    filt_feature = nrow(filt_mat),
    raw_bc = ncol(raw_mat),
    filt_bc = ncol(filt_mat)
  ))

  # Remove Raw Matrix to save space
  raw_mat <- NULL

  # Renaming the features
  filt_mat@Dimnames[[1]] <- str_replace(filt_mat@Dimnames[[1]], "_", "-")

  # Create Seurat Object
  sob.raw <- CreateSeuratObject(
    counts = filt_mat, min.cells = 200,
    min.features = 1000, project = i
  )

  cat(paste0("\nCheck-2 (", i, "): ", "Created Seurat Object\n"))

  # Check
  match_genes <- length(rownames(sob.raw)[(rownames(sob.raw) %in% mds_genes)])
  if (match_genes == length(mds_genes)) {
    cat(paste0("\nCheck-3 (", i, "): All Genes exist\n"))
  } else if (match_genes < length(mds_genes)) {
    rem <- length(mds_genes) - match_genes
    cat(paste0("\nCheck-3 (", i, "): ", rem, " Dropped\n"))
  }

  # Save the raw object
  # file_name <- paste0(prefixOut, i, "_raw_sob")
  # SaveH5Seurat(object = sob.raw, filename = file_name, overwrite = T,
  #              verbose = FALSE)

  # Calculate Percentage of Mitochondrial Read
  sob.raw[["percent.mt"]] <- PercentageFeatureSet(sob.raw, pattern = "^MT-")

  # Plots
  p <- VlnPlot(sob.raw,
    features = c("nFeature_RNA"),
    pt.size = 1
  ) +
    theme(
      title = element_text(size = 25),
      axis.text = element_text(size = 25),
      legend.position = "none"
    ) + xlab(paste0(
      "Individual: ", individual,
      " Age: ", age,
      " Sex: ", sex
    ))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Raw_VlnPlot_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- VlnPlot(sob.raw,
    features = c("nCount_RNA"),
    pt.size = 1
  ) +
    theme(
      title = element_text(size = 25),
      axis.text = element_text(size = 25),
      legend.position = "none"
    ) + xlab(paste0(
      "Individual: ", individual,
      " Age: ", age,
      " Sex: ", sex
    ))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Raw_VlnPlot_nCount_RNA.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- VlnPlot(sob.raw,
    features = c("percent.mt"),
    pt.size = 1
  ) +
    theme(
      title = element_text(size = 25),
      axis.text = element_text(size = 25),
      legend.position = "none"
    ) + xlab(paste0(
      "Individual: ", individual,
      " Age: ", age,
      " Sex: ", sex
    ))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Raw_VlnPlot_percent.mt.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme(
      title = element_text(size = 25), axis.text = element_text(size = 25),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Raw_FeatureScatter_percent.mt.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(
      title = element_text(size = 25), axis.text = element_text(size = 20),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Raw_FeatureScatter_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE
  )

  # Subset # Check Required if-else
  sob.sub <- subset(sob.raw, subset = nFeature_RNA > 200 & nCount_RNA < 30000 & percent.mt < 5)

  # Add Measuremnets
  measure <- append(measure, list(
    sub_feature = nrow(sob.sub),
    sub_bc = ncol(sob.sub)
  ))

  file_name <- paste0(prefixOut, i, "_sub_sob")
  SaveH5Seurat(
    object = sob.sub, filename = file_name, overwrite = T,
    verbose = FALSE
  )

  cat(paste0("\nCheck-4 (", i, "): Subset Written\n"))

  # Plots
  p <- VlnPlot(sob.sub, features = c("nFeature_RNA")) +
    theme(
      title = element_text(size = 25),
      axis.text = element_text(size = 25),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Sub_VlnPlot_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- VlnPlot(sob.sub, features = c("nCount_RNA")) +
    theme(
      title = element_text(size = 25),
      axis.text = element_text(size = 25),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Sub_VlnPlot_nCount_RNA.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- VlnPlot(sob.sub, features = c("percent.mt")) +
    theme(
      title = element_text(size = 25),
      axis.text = element_text(size = 25),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Sub_VlnPlot_percent.mt.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- FeatureScatter(sob.sub, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme(
      title = element_text(size = 25), axis.text = element_text(size = 20),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Sub_FeatureScatter_percent.mt.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- FeatureScatter(sob.sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(
      title = element_text(size = 25), axis.text = element_text(size = 20),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
  ggsave(p,
    filename = paste0(prefixOut, i, "_Sub_FeatureScatter_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE
  )

  # Check
  match_genes <- length(rownames(sob.sub)[(rownames(sob.sub) %in% mds_genes)])
  if (match_genes == length(mds_genes)) {
    cat(paste0("\nCheck-5 (", i, "): All Genes exist\n"))
  } else if (match_genes < length(mds_genes)) {
    rem <- length(mds_genes) - match_genes
    cat(paste0("\nCheck-5 (", i, "): ", rem, " Dropped\n"))
  }

  # Pre-processing
  sob.prs <- NormalizeData(sob.sub, verbose = F)
  sob.prs <- FindVariableFeatures(sob.prs, selection.method = "vst", verbose = F)
  sob.prs <- ScaleData(sob.prs, features = rownames(sob.prs), verbose = F)
  sob.prs <- RunPCA(sob.prs,
    features = VariableFeatures(sob.prs),
    verbose = F, ndims.print = 0, nfeatures.print = 0
  )

  # Basic UMAP
  p <- DimPlot(sob.prs, reduction = "pca", pt.size = 1) + theme(legend.position = "none") +
    ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
  ggsave(p,
    filename = paste0(prefixOut, i, "_basic_PCA.png"),
    dpi = 1200, limitsize = FALSE
  )

  # Save
  file_name <- paste0(prefixOut, i, "_prep_sob")
  SaveH5Seurat(
    object = sob.prs, filename = file_name, overwrite = T,
    verbose = FALSE
  )

  cat(paste0("\nCheck-6 (", i, "): Pre-processed Written\n"))

  # Check how many genes are present in the dataset
  indata <- rownames(sob.prs)[rownames(sob.prs) %in% reg.out]

  # Keep only those which are present in data
  biomart.anno <- biomart.anno[biomart.anno$SYMBOL %in% indata, ]

  # For seurat we need to divide genes into vectors
  s.genes <- unique(biomart.anno[biomart.anno$GO %in% c("GO:0006260"), "SYMBOL"])
  g2m.genes <- unique(biomart.anno[biomart.anno$GO %in% c("GO:0000087", "GO:0000279", "GO:0007059", "GO:0048285"), "SYMBOL"])

  # Cell Cycle Scores
  sob.cc <- CellCycleScoring(sob.prs,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )

  # Re-Run PCA
  sob.cc <- RunPCA(sob.cc,
    features = c(s.genes, g2m.genes),
    verbose = F, ndims.print = 0, nfeatures.print = 0
  )

  # Save The plot
  p <- DimPlot(sob.cc, pt.size = 1) +
    ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
    scale_color_manual(
      name = "Cell-Cycle GOs",
      breaks = c("G2M", "S", "G1"),
      values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
    )
  ggsave(p,
    filename = paste0(prefixOut, i, "_CCG_PCA.png"),
    dpi = 1200, limitsize = FALSE
  )

  # Scale to Regress out
  sob.cc <- ScaleData(sob.cc,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = rownames(sob.cc), verbose = F
  )

  # Re-Run PCA
  test <- RunPCA(sob.cc,
    features = c(s.genes, g2m.genes),
    verbose = F, ndims.print = 0, nfeatures.print = 0
  )
  p <- DimPlot(test, reduction = "pca", pt.size = 1) +
    ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
    scale_color_manual(
      name = "Cell-Cycle GOs",
      breaks = c("G2M", "S", "G1"),
      values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
    )
  ggsave(p,
    filename = paste0(prefixOut, i, "_CCC_PCA.png"),
    dpi = 1200, limitsize = FALSE
  )

  sob.cc <- RunPCA(sob.cc,
    features = VariableFeatures(sob.cc),
    nfeatures.print = 0, verbose = F, ndims.print = 0
  )

  p <- DimPlot(sob.cc, reduction = "pca", pt.size = 1) +
    ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
    scale_color_manual(
      name = "Cell-Cycle GOs",
      breaks = c("G2M", "S", "G1"),
      values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
    )

  ggsave(p,
    filename = paste0(prefixOut, i, "_CCC_all_PCA.png"),
    dpi = 1200, limitsize = FALSE
  )

  # Save
  file_name <- paste0(prefixOut, i, "_CCC_sob")
  SaveH5Seurat(
    object = sob.cc, filename = file_name, overwrite = T,
    verbose = FALSE
  )

  cat(paste0("\nCheck-7 (", i, "): Cell-Cycle-Corrected\n"))

  # More Pre-processing
  sob.p <- FindNeighbors(sob.cc, verbose = F)
  sob.p <- FindClusters(sob.p, verbose = F)
  sob.p <- RunUMAP(sob.p, dims = 1:50, verbose = F)

  cat(paste0("\nCheck-8 (", i, "): Clustering Completed\n"))

  file_name <- paste0(prefixOut, i, "_Processed_sob")
  SaveH5Seurat(
    object = sob.p, filename = file_name,
    overwrite = T, verbose = FALSE
  )

  # Save the UMAP with the cluster
  p <- DimPlot(sob.p, reduction = "umap", pt.size = 1, group.by = "seurat_clusters") +
    ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
    scale_color_hue(l = 50)
  ggsave(p,
    filename = paste0(prefixOut, i, "_C_UMAP.png"),
    dpi = 1200, limitsize = FALSE
  )

  cat(paste0("\nCheck-9 (", i, "): Saved Processed File\n"))

  # Save
  write.table(t(as.data.frame(measure)),
    paste0(prefixOut, i, "_Measure.tsv"),
    col.names = NA
  )

  cat(paste0("\nCheck-10: Completed for ", i, "\n"))
}
