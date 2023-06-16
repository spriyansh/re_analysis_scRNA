# Call the required libraries
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))

# Trend Plot themes
trend_theme <- theme(
  plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
  legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
  axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
  legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
)

umap_theme <- theme(
  plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
  legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
  axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
  legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
)

# Prefix
prefixIn <- "/home/priyansh/gitDockers/re_analysis_scRNA/setty_et_al/data/output/monocle3/Monocle3UMAPOB/"
dir.create("/home/priyansh/gitDockers/re_analysis_scRNA/setty_et_al/data/output/monocle3/MarkersTrends", showWarnings = FALSE)
prefixOut <- "/home/priyansh/gitDockers/re_analysis_scRNA/setty_et_al/data/output/monocle3/MarkersTrends/"
prefixOut2 <- "/home/priyansh/gitDockers/re_analysis_scRNA/setty_et_al/data/output/monocle3/MarkerUmaps/"

# mds genes
mds_genes <- c(
  "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
  "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
  "NF1", "CALR", "MPL", "GATA2"
)

# Other Explorative Genes
oth_genes <- c("MPO", "CD34", "CD79B", "GATA1", "IRF8", "CD41", "EPOR", "EPOR", "AICDA", "APOBEC3C", "APOBEC3F", "APOBEC3G")

# all_file_names
all_file_names <- list.files(paste0(prefixIn))
names(all_file_names) <- str_remove(all_file_names, "_m3_UMAP_OB.RDS")

# Running for each object
for (i in names(all_file_names)) {
  if (i == "rep1") {
    age <- "35"
    sex <- "Male"
  } else if (i == "rep2") {
    age <- "28"
    sex <- "Female"
  } else if (i == "rep3") {
    age <- "19"
    sex <- "Female"
  }

  # Get filename
  fileName <- all_file_names[[i]]

  # Load object
  cds <- readRDS(paste0(prefixIn, fileName))

  # Subset the marker genes
  marker.genes <- c(mds_genes, oth_genes)

  # Subset aviable markers
  marker.genes <- marker.genes[marker.genes %in% rowData(cds)$gene_short_name]

  for (j in c(marker.genes)) {
    trend.plot <- plot_genes_in_pseudotime(
      cds_subset = cds[rowData(cds)$gene_short_name == j, ],
      min_expr = 0, color_cells_by = "cell.type",
      label_by_short_name = F
    ) + scale_color_brewer(palette = "Dark2") +
      ggtitle(paste("Marker Genes: ", j),
        subtitle = paste("Age", age, "Sex", sex)
      )
    trend.plot <- trend.plot + trend_theme
    trend.plot <- trend.plot + guides(colour = guide_legend(override.aes = list(size = rel(6))))
    ggsave(trend.plot,
      filename = paste0(prefixOut, i, "_", j, "_Marker.png"),
      dpi = 1300, limitsize = FALSE, width = 6
    )
  }

  for (j in c(marker.genes)) {
    marker.umap <- plot_cells(cds,
      reduction_method = "UMAP", cell_size = 1,
      label_cell_groups = F, show_trajectory_graph = F,
      genes = j
    ) + viridis::scale_color_viridis(option = "E", direction = -1) +
      ggtitle(paste("Lineage Marker Gene: ", j),
        subtitle = paste("Age", age, "Sex", sex)
      )
    marker.umap <- marker.umap + umap_theme
    ggsave(marker.umap,
      filename = paste0(prefixOut2, i, "_", j, ".png"),
      dpi = 1300, limitsize = FALSE
    )
  }
}
