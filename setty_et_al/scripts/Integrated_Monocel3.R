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
prefixIn <- "setty_et_al/data/output/Integrated/"
prefixOut <- "setty_et_al/data/output/Integrated/"
sob <- LoadH5Seurat(file = paste0(prefixIn, "Setty_Integrated_sob.h5seurat"), verbose = F)

# CDS Temp
cds_temp <- new_cell_data_set(
  expression_data = sob@assays$RNA@data,
  cell_metadata = sob@meta.data,
  gene_metadata = data.frame(
    gene_short_name = rownames(sob@assays$RNA@data),
    row.names = rownames(sob@assays$RNA@data)
  )
)
cds_temp <- preprocess_cds(cds_temp,
  num_dim = 10, method = "PCA",
  norm_method = "none"
)
cds_temp <- reduce_dimension(cds_temp,
  preprocess_method = "PCA",
  reduction_method = "UMAP", umap.min_dist = 0.1,
  umap.n_neighbors = 20
)

# Plotting Via Donors
p <- plot_cells(cds_temp,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "orig.ident", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
  ggtitle("All Individual Together") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Unintegrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Unintegrated/Donor.png"),
  dpi = 1300, limitsize = FALSE
)

# Plotting Via Sex
p <- plot_cells(cds_temp,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "sex", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
  ggtitle("All Individual Together") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Unintegrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Unintegrated/Sex.png"),
  dpi = 1300, limitsize = FALSE
)

# Plotting Via Age
p <- plot_cells(cds_temp,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "age", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
  ggtitle("All Individual Together") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Unintegrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Unintegrated/Age.png"),
  dpi = 1300, limitsize = FALSE
)

# Plotting Via Cell.Type
p <- plot_cells(cds_temp,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "cell.type", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
  ggtitle("All Individual Together") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Unintegrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Unintegrated/CellType.png"),
  dpi = 1300, limitsize = FALSE
)

# CDS Integrated
cds <- new_cell_data_set(
  expression_data = sob@assays$integrated@data,
  cell_metadata = sob@meta.data,
  gene_metadata = data.frame(
    gene_short_name = rownames(sob@assays$integrated@data),
    row.names = rownames(sob@assays$integrated@data)
  )
)
cds <- preprocess_cds(cds,
  num_dim = 10, method = "PCA",
  norm_method = "none"
)
cds <- reduce_dimension(cds,
  preprocess_method = "PCA",
  reduction_method = "UMAP", umap.min_dist = 0.3,
  umap.n_neighbors = 25
)

# Plotting Via Donors
p <- plot_cells(cds,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "orig.ident", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
  ggtitle("Integrated Data") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Integrated/Donor.png"),
  dpi = 1300, limitsize = FALSE
)

# Plotting Via Sex
p <- plot_cells(cds,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "sex", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
  ggtitle("Integrated Data") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Integrated/Sex.png"),
  dpi = 1300, limitsize = FALSE
)

# Plotting Via Age
p <- plot_cells(cds,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "age", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
  ggtitle("Integrated Data") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Integrated/Age.png"),
  dpi = 1300, limitsize = FALSE
)

# Plotting Via Cell.Type
p <- plot_cells(cds,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "cell.type", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
  ggtitle("Integrated Data") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Integrated/CellType.png"),
  dpi = 1300, limitsize = FALSE
)

# Perform Liden Clustering
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Plotting Monocle3 UMAP
p <- plot_cells(cds,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "cluster", label_cell_groups = F
) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Integrated Data") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p,
  filename = paste0(prefixOut, "/Integrated/Cluster.png"),
  dpi = 1300, limitsize = FALSE
)

# Learn Graph
cds <- learn_graph(cds, verbose = F)

# Plotting Monocle3 UMAP
p1 <- plot_cells(cds,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "cell.type", label_cell_groups = F
) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Integrated Data") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  ) +
  guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p1,
  filename = paste0(prefixOut, "/Integrated/TrLine.png"),
  dpi = 1300, limitsize = FALSE
)

# Order cells
cds <- order_cells(cds, root_cells = rownames(
  colData(cds)[(colData(cds)$cell.type == "HSC" &
    colData(cds)$predicted.celltype.l2.score > 0.9 &
    colData(cds)$mapping.score > 0.9), ]
))

# Plotting Monocle3 UMAP
p2 <- plot_cells(cds,
  reduction_method = "UMAP", cell_size = 0.5,
  color_cells_by = "pseudotime", label_cell_groups = F,
  label_branch_points = F, label_leaves = F,
  label_roots = T, label_principal_points = F
) +
  viridis::scale_color_viridis(option = "C") +
  ggtitle("Integrated Data") +
  theme(
    plot.margin = margin(0.2, 0.2, 0, 0, "cm"),
    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))
  )
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p2,
  filename = paste0(prefixOut, "/Integrated/Pseudotime.png"),
  dpi = 1300, limitsize = FALSE
)

p <- p1 + p2
ggsave(p,
  filename = paste0(prefixOut, "/Integrated/PseudotimeCombined.png"),
  dpi = 1300, limitsize = FALSE, width = 8
)

# Save the Object
saveRDS(cds, paste0(prefixOut, "/Integrated/m3.RDS"))
