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
prefixIn <- "guess_et_al/data/input/p3/"
prefixOut <- "guess_et_al/data/output/p3/"

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
cite_file_names <- all_file_names[all_file_names != grep(all_file_names, pattern = "h5", value = T)] 

# Seurat H5
cr_counts <- grep(all_file_names, pattern = "h5", value = T) 
cr_counts <- Read10X_h5(paste0(prefixIn, cr_counts))

# Read all the HTO
cite.list <- list()

# Read
for (i in cite_file_names){
    matrix_dir = paste0(prefixIn, i)
    barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
    features.path <- paste0(matrix_dir, "/features.tsv.gz")
    matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = feature.names$V1
    sob <- CreateSeuratObject(counts = mat, project = i)
    cite.list <- append(cite.list, list(sob))
}

# ADD names
names(cite.list) <- cite_file_names

# Intersect Barcodes
joint.bcs <- intersect(colnames(cr_counts),
                       colnames(cite.list[[1]]),
                       colnames(cite.list[[2]]),
                       colnames(cite.list[[3]]),
                       colnames(cite.list[[4]]),
                       colnames(cite.list[[5]]),
                       colnames(cite.list[[6]]),
                       colnames(cite.list[[7]]),
                       colnames(cite.list[[8]]))


# Distribute filenames
MDS17_names <- grep(all_file_names, pattern = "MDS", value = T) # Male-35
AML17_names <- grep(all_file_names, pattern = "AML", value = T) # Female-28

# Make a list of filenames
rep_list <- list(mds17 = MDS17_names, aml17 = AML17_names)

# Run for every dataset
for (i in names(rep_list)) {
  # i <- "aml17"

  if (i == "mds17") {
    patientID <- "17"
    condition <- "MDS"
  } else if (i == "aml17") {
      patientID <- "17"
      condition <- "sAML"
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
    counts = filt_mat, min.cells = 500,
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
    pt.size = 1
  ) +
    theme(
      title = element_text(size = 25),
      axis.text = element_text(size = 25),
      legend.position = "none"
    ) + xlab(paste0(
      "PatientID: ", patientID,
      " Condition: ", condition
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_VlnPlot_nFeature_RNA.png"),
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
        "PatientID: ", patientID,
        " Condition: ", condition
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_VlnPlot_nCount_RNA.png"),
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
        "PatientID: ", patientID,
        " Condition: ", condition
    ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_VlnPlot_percent.mt.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    theme(
      title = element_text(size = 25), axis.text = element_text(size = 25),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
      xlab(paste0("PatientID: ", patientID," Condition: ", condition))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_FeatureScatter_percent.mt.png"),
    dpi = 1200, limitsize = FALSE
  )

  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(
      title = element_text(size = 25), axis.text = element_text(size = 20),
      legend.position = "none"
    ) + geom_smooth(method = "lm", formula = "y~x") +
      xlab(paste0("PatientID: ", patientID," Condition: ", condition))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  ggsave(p,
    filename = paste0(outDir, i, "_Raw_FeatureScatter_nFeature_RNA.png"),
    dpi = 1200, limitsize = FALSE
  )

  # Subset # Check Required if-else
  sob.sub <- subset(sob.raw, subset = nFeature_RNA > 100 & nCount_RNA < 30000 & percent.mt < 15)

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
               pt.size = 1
  ) +
      theme(
          title = element_text(size = 25),
          axis.text = element_text(size = 25),
          legend.position = "none"
      ) + xlab(paste0(
          "PatientID: ", patientID,
          " Condition: ", condition
      ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
         filename = paste0(outDir, i, "_Sub_VlnPlot_nFeature_RNA.png"),
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
          "PatientID: ", patientID,
          " Condition: ", condition
      ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
         filename = paste0(outDir, i, "_Sub_VlnPlot_nCount_RNA.png"),
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
          "PatientID: ", patientID,
          " Condition: ", condition
      ))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
         filename = paste0(outDir, i, "_Sub_VlnPlot_percent.mt.png"),
         dpi = 1200, limitsize = FALSE
  )
  
  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      theme(
          title = element_text(size = 25), axis.text = element_text(size = 25),
          legend.position = "none"
      ) + geom_smooth(method = "lm", formula = "y~x") +
      xlab(paste0("PatientID: ", patientID," Condition: ", condition))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  ggsave(p,
         filename = paste0(outDir, i, "_Sub_FeatureScatter_percent.mt.png"),
         dpi = 1200, limitsize = FALSE
  )
  
  p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      theme(
          title = element_text(size = 25), axis.text = element_text(size = 20),
          legend.position = "none"
      ) + geom_smooth(method = "lm", formula = "y~x") +
      xlab(paste0("PatientID: ", patientID," Condition: ", condition))
  outDir <- paste0(prefixOut, "p17/QC_Plots/")
  ggsave(p,
         filename = paste0(outDir, i, "_Sub_FeatureScatter_nFeature_RNA.png"),
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
  sob.prs <- FindVariableFeatures(sob.prs, selection.method = "vst", verbose = F,
                                  nfeatures = 6000)
  sob.prs <- ScaleData(sob.prs, features = rownames(sob.prs), verbose = F)
  sob.prs <- RunPCA(sob.prs,
    features = VariableFeatures(sob.prs),
    verbose = F, ndims.print = 0, nfeatures.print = 0
  )

  # Basic UMAP
  p <- DimPlot(sob.prs, reduction = "pca", pt.size = 1) + theme(legend.position = "none") +
    ggtitle(paste0("PatientID: ", patientID," Condition: ", condition))
  outDir <- paste0(prefixOut, "p17/PCA/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_basic_PCA.png"),
    dpi = 1200, limitsize = FALSE
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
  sob.cc <- RunPCA(sob.cc,npcs = 10,
    features = c(s.genes, g2m.genes),
    verbose = F, ndims.print = 0, nfeatures.print = 0
  )

  # Save The plot
  p <- DimPlot(sob.cc, pt.size = 1) +
    ggtitle(paste0("PatientID: ", patientID," Condition: ", condition)) +
    scale_color_manual(
      name = "Cell-Cycle GOs",
      breaks = c("G2M", "S", "G1"),
      values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
    )
  outDir <- paste0(prefixOut, "p17/PCA/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
    filename = paste0(outDir, i, "_CCG_PCA.png"),
    dpi = 1200, limitsize = FALSE
  )

  # Scale to Regress out
  sob.cc <- ScaleData(sob.cc,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = rownames(sob.cc), verbose = T
  )

  # Re-Run PCA
  test <- RunPCA(sob.cc,npcs = 10,
    features = c(s.genes, g2m.genes),
    verbose = F, ndims.print = 0, nfeatures.print = 0
  )
  p <- DimPlot(test, reduction = "pca", pt.size = 1) +
      ggtitle(paste0("PatientID: ", patientID," Condition: ", condition)) +
    scale_color_manual(
      name = "Cell-Cycle GOs",
      breaks = c("G2M", "S", "G1"),
      values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
    )
  outDir <- paste0(prefixOut, "p17/PCA/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
         filename = paste0(outDir, i, "_CCC_PCA.png"),
         dpi = 1200, limitsize = FALSE
  )

  sob.cc <- RunPCA(sob.cc,
    features = VariableFeatures(sob.cc),
    nfeatures.print = 0, verbose = F, ndims.print = 0
  )

  p <- DimPlot(sob.cc, reduction = "pca", pt.size = 1) +
      ggtitle(paste0("PatientID: ", patientID," Condition: ", condition)) +
    scale_color_manual(
      name = "Cell-Cycle GOs",
      breaks = c("G2M", "S", "G1"),
      values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
    )
  outDir <- paste0(prefixOut, "p17/PCA/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
         filename = paste0(outDir, i, "_CCC_all_PCA.png"),
         dpi = 1200, limitsize = FALSE
  )

  # Save
  outDir <- paste0(prefixOut, "p17/SeuratObjects/")
  dir.create(outDir, showWarnings = F)
  file_name <- paste0(outDir, i, "_CCC_sob")
  SaveH5Seurat(
      object = sob.prs, filename = file_name, overwrite = T,
      verbose = FALSE
  )

  cat(paste0("\nCheck-7 (", i, "): Cell-Cycle-Corrected\n"))

  # More Pre-processing
  sob.p <- FindNeighbors(sob.cc, verbose = F)
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
  p <- DimPlot(sob.p, reduction = "umap", pt.size = 1, group.by = "seurat_clusters") +
      ggtitle(paste0("PatientID: ", patientID," Condition: ", condition)) +
    scale_color_hue(l = 50)
  outDir <- paste0(prefixOut, "p17/UMAP/")
  dir.create(outDir, showWarnings = F)
  ggsave(p,
         filename = paste0(outDir, i, "_C_UMAP.png"),
         dpi = 1200, limitsize = FALSE
  )

  cat(paste0("\nCheck-9 (", i, "): Saved Processed File\n"))

  # Save
  write.table(t(as.data.frame(measure)),
    paste0(prefixOut, "/p17/",i, "p17_Measure.tsv"),
    col.names = NA
  )

  cat(paste0("\nCheck-10: Completed for ", i, "\n"))
}
