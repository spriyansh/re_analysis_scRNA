#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Monocle3 Integrate #########
#######################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))

# Prefix
prefixIn <- "setty_et_al/data/output/Integrated/Integrated/"
prefixOut <- "setty_et_al/data/output/Integrated/"
cds <- readRDS(file = paste0(prefixIn, "m3.RDS"))

# Run Autocorrelation
auto.res <- graph_test(cds, neighbor_graph = "principal_graph")

# Save
write.table(as.data.frame(auto.res),
  sep = "\t",
  quote = F, row.names = F, file = paste0(prefixOut, "/AutoCorr.tsv")
)

# Transitory Genes
trans_genes <- auto.res[auto.res$status == "OK", "gene_short_name"]

# mds genes
mds_genes <- c(
  "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
  "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
  "NF1", "CALR", "MPL", "GATA2"
)

# Other Explorative Genes
LM_genes <- c("MPO", "CD34", "CD79B", "GATA1", "IRF8", "CD41", "EPOR", "EPOR", "AICDA", "APOBEC3C", "APOBEC3F", "APOBEC3G")

# Selection
mds_genes <- mds_genes[mds_genes %in% rownames(cds)]
LM_genes <- LM_genes[LM_genes %in% rownames(cds)]

# Create Directories
dir.create(paste0(prefixOut, "/Transit_UMAP/"), showWarnings = FALSE)
dir.create(paste0(prefixOut, "/MDS_UMAP/"), showWarnings = FALSE)
dir.create(paste0(prefixOut, "/LM_UMAP/"), showWarnings = FALSE)
dir.create(paste0(prefixOut, "/Transit_Trend/"), showWarnings = FALSE)
dir.create(paste0(prefixOut, "/MDS_Trend/"), showWarnings = FALSE)
dir.create(paste0(prefixOut, "/LM_Trend/"), showWarnings = FALSE)

# Create Plots For Transitory Genes in UMAP
for (i in trans_genes) {
  # Plotting Via Donors
  p <- plot_cells(cds,
    reduction_method = "UMAP", cell_size = 0.5,
    label_cell_groups = F, show_trajectory_graph = F,
    genes = i
  ) +
    viridis::scale_color_viridis(option = "E", direction = -1) +
    ggtitle(paste("Transitory Gene: ", i)) +
    theme(
      plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
      legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
      axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
      legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
    )
  ggsave(p,
    filename = paste0(prefixOut, "/Transit_UMAP/", i, ".png"),
    dpi = 1300, limitsize = FALSE
  )
}

# Create Plots For Transitory Genes in UMAP
for (i in mds_genes) {
  # Plotting Via Donors
  p <- plot_cells(cds,
    reduction_method = "UMAP", cell_size = 0.5,
    label_cell_groups = F, show_trajectory_graph = F,
    genes = i
  ) +
    viridis::scale_color_viridis(option = "E", direction = -1) +
    ggtitle(paste("MDS Gene: ", i)) +
    theme(
      plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
      legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
      axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
      legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
    )
  ggsave(p,
    filename = paste0(prefixOut, "/MDS_UMAP/", i, ".png"),
    dpi = 1300, limitsize = FALSE
  )
}


# Create Plots For Transitory Genes in UMAP
for (i in LM_genes) {
  # Plotting Via Donors
  p <- plot_cells(cds,
    reduction_method = "UMAP", cell_size = 0.5,
    label_cell_groups = F, show_trajectory_graph = F,
    genes = i
  ) +
    viridis::scale_color_viridis(option = "E", direction = -1) +
    ggtitle(paste("Lineage Marker Gene: ", i)) +
    theme(
      plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
      legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
      axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
      legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
    )
  ggsave(p,
    filename = paste0(prefixOut, "/LM_UMAP/", i, ".png"),
    dpi = 1300, limitsize = FALSE
  )
}

# Create Plots For Transitory Genes in UMAP
for (i in trans_genes) {
  # Plotting Via Donors
  p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == i, ],
    min_expr = 0.5, color_cells_by = "cell.type"
  ) +
    scale_color_brewer(palette = "Dark2") +
    ggtitle(paste("Transitory Gene: ", i)) +
    theme(
      plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
      legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
      axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
      legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
    ) +
    guides(colour = guide_legend(override.aes = list(size = rel(6))))
  ggsave(p,
    filename = paste0(prefixOut, "/Transit_Trend/", i, ".png"),
    dpi = 1300, limitsize = FALSE
  )
}

# Create Plots For Transitory Genes in UMAP
for (i in mds_genes) {
  # Plotting Via Donors
  p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == i, ],
    min_expr = 0.5, color_cells_by = "cell.type"
  ) +
    scale_color_brewer(palette = "Dark2") +
    ggtitle(paste("MDS Gene: ", i)) +
    theme(
      plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
      legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
      axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
      legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
    ) +
    guides(colour = guide_legend(override.aes = list(size = rel(6))))
  ggsave(p,
    filename = paste0(prefixOut, "/MDS_Trend/", i, ".png"),
    dpi = 1300, limitsize = FALSE
  )
}

# Create Plots For Transitory Genes in UMAP
for (i in LM_genes) {
  # Plotting Via Donors
  p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == i, ],
    min_expr = 0, color_cells_by = "cell.type",
    label_by_short_name = F
  ) +
    scale_color_brewer(palette = "Dark2") +
    ggtitle(paste("Lineage Marker Gene: ", i)) +
    theme(
      plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
      legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
      axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
      legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
    ) +
    guides(colour = guide_legend(override.aes = list(size = rel(6))))
  ggsave(p,
    filename = paste0(prefixOut, "/LM_Trend/", i, ".png"),
    dpi = 1300, limitsize = FALSE
  )
}

#
# spline_calc <- na.omit(data.frame(pseudotime = pseudotime(cds),
#                                   counts =  as.vector(assay(cds_sub)),
#                                   cell.type = colData(cds_sub)$cell.type))
#
# spline_calc <- spline_calc[!is.na(spline_calc$pseudotime), ]
# spline_calc <- spline_calc[spline_calc$pseudotime != Inf, ]
#
# splineDf <- ns(x = spline_calc$counts, df = 3, knots = 2)
#
# splineDf <- cbind(splineDf, spline_calc)
# colnames(splineDf) <- c("x_spline", "y_spline", "pseudotime", "count", "cell.type")
#
# # Plot With ggplot2
# ggplot(splineDf, aes(x = pseudotime, y = count, color = cell.type))+
#     geom_point(size = 1, alpha = 0.8) + scale_color_brewer(palette = "Dark2") +
#     geom_line(aes(x= x_spline, y = y_spline), color = "black", linewidth = 3)+
#     ggtitle(paste("Lineage Marker Gene: ", i)) + theme_classic()+
#     theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
#           legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
#           axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
#           legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
#     guides(colour = guide_legend(override.aes = list(size = rel(6))))
