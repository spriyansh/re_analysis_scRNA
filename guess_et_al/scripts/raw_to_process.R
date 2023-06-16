#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Pre-Process ################
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
prefixIn <- "guess_et_al/data/input/"
prefixOut <- "guess_et_al/data/output/"

# mds genes
mds_genes <- c(
  "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
  "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
  "NF1", "CALR", "MPL", "GATA2"
)


# all_file_names
all_file_names <- list.files(paste0(prefixIn, "/p17/raw_barcode_matrix/"))

# Distribute filenames
MDS17_names <- grep(all_file_names, pattern = "MDS", value = T) # Male-84-MDS
AML17_names <- grep(all_file_names, pattern = "AML", value = T) # Male-84-AML

# Make a list of filenames
rep_list <- list(MDS17 = MDS17_names, AML17 = AML17_names)

# Run for every dataset
for (i in names(rep_list)) {
  # i <- "MDS17"

  if (i == "MDS17") {
    patientID <- "17"
    condition <- "MDS"
    age <- 84
    sex <- "Male"
  } else if (i == "AML17") {
    patientID <- "17"
    condition <- "sAML"
    age <- 84
    sex <- "Male"
  }

  measure <- list()

  # Get the names rep#
  rep_i_name <- rep_list[[i]]

  # Load the RDS files
  filt_mat <- Read10X_h5(paste0(prefixIn, "/p17/raw_barcode_matrix/", rep_i_name))

  cat(paste0("\nCheck-1 (", i, "): ", "Loaded Files\n"))

  # Add Measures to the list
  measure <- append(measure, list(
    filt_feature = nrow(filt_mat),
    filt_bc = ncol(filt_mat)
  ))

  # Renaming the features
  filt_mat@Dimnames[[1]] <- str_replace(filt_mat@Dimnames[[1]], "_", "-")

  # Create Seurat Object
  sob.raw <- CreateSeuratObject(
    counts = filt_mat, min.cells = 250,
    min.features = 250, project = i
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
    pt.size = 0.5
  ) +
    theme(
      title = element_text(size = 10),
      axis.text = element_text(size = 7),
      legend.position = "none"
    ) + xlab(paste0(
      "PatientID: ", patientID,
      " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_VlnPlot_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  p <- VlnPlot(sob.raw,
    features = c("nCount_RNA"),
    pt.size = 0.5
  ) +
    theme(
      title = element_text(size = 10),
      axis.text = element_text(size = 7),
      legend.position = "none"
    ) + xlab(paste0(
      "PatientID: ", patientID,
      " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_VlnPlot_nCount_RNA.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  p <- VlnPlot(sob.raw,
    features = c("percent.mt"),
    pt.size = 0.5
  ) +
    theme(
      title = element_text(size = 10),
      axis.text = element_text(size = 7),
      legend.position = "none"
    ) + xlab(paste0(
      "PatientID: ", patientID,
      " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_VlnPlot_percent.mt.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme(
      title = element_text(size = 10), axis.text = element_text(size = 7),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0(
      "PatientID: ", patientID, " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_FeatureScatter_percent.mt.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(
      title = element_text(size = 10), axis.text = element_text(size = 7),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0(
      "PatientID: ", patientID, " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_FeatureScatter_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  # Subset # Check Required if-else
  sob.sub <- subset(sob.raw, subset = percent.mt < 15)

  # Add Measuremnets
  measure <- append(measure, list(
    sub_feature = nrow(sob.sub),
    sub_bc = ncol(sob.sub)
  ))

  outDir <- paste0(prefixOut, "p17/SeuratObjects/")
  dir.create(outDir, showWarnings = F)
  file_name <- paste0(outDir, i, "_sub_sob")
  SaveH5Seurat(
    object = sob.sub, filename = file_name, overwrite = T,
    verbose = FALSE
  )

  cat(paste0("\nCheck-4 (", i, "): Subset Written\n"))

  # Plots
  p <- VlnPlot(sob.raw,
    features = c("nFeature_RNA"),
    pt.size = 0.5
  ) +
    theme(
      title = element_text(size = 10),
      axis.text = element_text(size = 7),
      legend.position = "none"
    ) + xlab(paste0(
      "PatientID: ", patientID,
      " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Sub_VlnPlot_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  p <- VlnPlot(sob.raw,
    features = c("nCount_RNA"),
    pt.size = 0.5
  ) +
    theme(
      title = element_text(size = 10),
      axis.text = element_text(size = 7),
      legend.position = "none"
    ) + xlab(paste0(
      "PatientID: ", patientID,
      " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Sub_VlnPlot_nCount_RNA.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  p <- VlnPlot(sob.raw,
    features = c("percent.mt"),
    pt.size = 0.5
  ) +
    theme(
      title = element_text(size = 10),
      axis.text = element_text(size = 7),
      legend.position = "none"
    ) + xlab(paste0(
      "PatientID: ", patientID,
      " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Sub_VlnPlot_percent.mt.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme(
      title = element_text(size = 10), axis.text = element_text(size = 7),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0(
      "PatientID: ", patientID, " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  ggsave(p,
    filename = paste0(outDir, i, "_Sub_FeatureScatter_percent.mt.png"),
    dpi = 1200, limitsize = FALSE, width = 4
  )

  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(
      title = element_text(size = 10), axis.text = element_text(size = 7),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0(
      "PatientID: ", patientID, " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  ggsave(p,
    filename = paste0(outDir, i, "_Sub_FeatureScatter_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE, width = 4
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
  sob.prs <- FindVariableFeatures(sob.prs,
    selection.method = "vst", verbose = F,
    nfeatures = 6000
  )
  sob.prs <- ScaleData(sob.prs, features = rownames(sob.prs), verbose = F)
  sob.prs <- RunPCA(sob.prs,
    features = VariableFeatures(sob.prs),
    verbose = F, ndims.print = 0, nfeatures.print = 0
  )

  # Basic UMAP
  p <- DimPlot(sob.prs, reduction = "pca", pt.size = 0.5) + theme(legend.position = "none") +
    ggtitle(paste0(
      "PatientID: ", patientID, " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    ))
  outDir <- paste0(prefixOut, "p17/PCA/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_basic_PCA.png"),
    dpi = 1200, limitsize = FALSE, width = 8
  )

  # Save
  outDir <- paste0(prefixOut, "p17/SeuratObjects/")
  dir.create(outDir, showWarnings = F)
  file_name <- paste0(outDir, i, "_prep_sob")
  SaveH5Seurat(
    object = sob.prs, filename = file_name, overwrite = T,
    verbose = FALSE
  )

  cat(paste0("\nCheck-6 (", i, "): Pre-processed Written\n"))

  # More Pre-processing
  sob.p <- FindNeighbors(sob.prs, verbose = F)
  sob.p <- FindClusters(sob.p, verbose = F)
  sob.p <- RunUMAP(sob.p, dims = 1:10, verbose = F)

  cat(paste0("\nCheck-8 (", i, "): Clustering Completed\n"))

  outDir <- paste0(prefixOut, "p17/SeuratObjects/")
  dir.create(outDir, showWarnings = F)
  file_name <- paste0(outDir, i, "_Processed_sob")
  SaveH5Seurat(
    object = sob.p, filename = file_name, overwrite = T,
    verbose = FALSE
  )

  # Save the UMAP with the cluster
  p <- DimPlot(sob.p, reduction = "umap", pt.size = 0.5, group.by = "seurat_clusters") +
    ggtitle(paste0(
      "PatientID: ", patientID, " Condition: ", condition,
      " Age: ", age,
      " Sex: ", sex
    )) +
    scale_color_hue(l = 50)
  outDir <- paste0(prefixOut, "p17/UMAP/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_C_UMAP.png"),
    dpi = 1200, limitsize = FALSE, width = 8
  )

  cat(paste0("\nCheck-9 (", i, "): Saved Processed File\n"))

  # Save
  write.table(t(as.data.frame(measure)),
    paste0(prefixOut, "/p17/", i, "p17_Measure.tsv"),
    col.names = NA
  )

  cat(paste0("\nCheck-10: Completed for ", i, "\n"))
}
