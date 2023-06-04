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
prefixIn <- "guess_et_al/data/output/p17/Integrated/"
prefixOut <- "guess_et_al/data/output/p17/Integrated/"
sob <- LoadH5Seurat(file = paste0(prefixIn, "Guess_Integrated_sob.h5seurat"), verbose = F)

# CDS Temp
cds_temp <- new_cell_data_set(
    expression_data = sob@assays$RNA@data,
    cell_metadata = sob@meta.data,
    gene_metadata = data.frame(
        gene_short_name = rownames(sob@assays$RNA@data),
        row.names = rownames(sob@assays$RNA@data)
    )
)
cds_temp <- preprocess_cds(cds_temp, num_dim = 10, method = "PCA",
                           norm_method = "none")
cds_temp <- reduce_dimension(cds_temp, preprocess_method = "PCA",
                             reduction_method = "UMAP", umap.min_dist = 0.1,
                             umap.n_neighbors = 20)

# Plotting Via Donors
p <- plot_cells(cds_temp, reduction_method = "UMAP", cell_size = 0.5,
                color_cells_by = "orig.ident", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
    ggtitle("Patient-17: MDS and sAML") + 
    theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
          legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
          axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
    guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Unintegrated/"), showWarnings = FALSE)
ggsave(p, filename = paste0(prefixOut, "/Unintegrated/Year.png"),
       dpi = 1300, limitsize = FALSE
)

# Plotting Via Cell.Type
p <- plot_cells(cds_temp, reduction_method = "UMAP", cell_size = 0.5,
                color_cells_by = "cell.type", label_cell_groups = F
) + #scale_color_brewer(palette = "Set1") +
    ggtitle("Patient-17: MDS and sAML") + 
    theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
          legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
          axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
    guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Unintegrated/"), showWarnings = FALSE)
ggsave(p, filename = paste0(prefixOut, "/Unintegrated/CellType.png"),
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
cds <- preprocess_cds(cds, num_dim = 10, method = "PCA",
                           norm_method = "none")
cds <- reduce_dimension(cds, preprocess_method = "PCA",
                             reduction_method = "UMAP", umap.min_dist = 0.3,
                             umap.n_neighbors = 25)

# Plotting Via Donors
p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 0.5,
                color_cells_by = "orig.ident", label_cell_groups = F
) + scale_color_brewer(palette = "Dark2") +
    ggtitle("Patient-17: MDS and sAML, Integrated Data") + 
    theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
          legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
          axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
    guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p, filename = paste0(prefixOut, "/Integrated/Donor.png"),
       dpi = 1300, limitsize = FALSE
)

# Plotting Via Cell.Type
p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 0.5,
                color_cells_by = "cell.type", label_cell_groups = F,
) + #scale_color_brewer(palette = "Dark2") +
    ggtitle("Patient-17: MDS and sAML, Integrated Data") + 
    theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
          legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
          axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
    guides(colour = guide_legend(override.aes = list(size = rel(6))))
dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p, filename = paste0(prefixOut, "/Integrated/CellType.png"),
       dpi = 1300, limitsize = FALSE
)


# Comparison
umap.df <- data.frame(cell.type = as.factor(colData(cds)$cell.type),
                      orig.ident = as.factor(colData(cds)$orig.ident),
                      umap1 = cds@reduce_dim_aux$UMAP$model$umap_model$embedding[,1],
                      umap2 = cds@reduce_dim_aux$UMAP$model$umap_model$embedding[,2])

p1 <- ggplot(umap.df, aes(x = umap1, y = umap2, alpha = orig.ident, 
                    color = cell.type)) +
    geom_point(data = umap.df, size = 0.5, shape = 16) + theme_classic() + ggtitle("Patient-17: MDS State") + 
    scale_alpha_manual(values = c(aml17 = 0.1, mds17 = 0.8))+
     theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
          legend.position = "none", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
          axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
    guides(colour = guide_legend(override.aes = list(size = rel(2))))

p2 <- ggplot(umap.df, aes(x = umap1, y = umap2, alpha = orig.ident, 
                          color = cell.type)) +
    geom_point(data = umap.df, size = 0.5, shape = 16) + theme_classic() + ggtitle("Patient-17: sAML State (After Two Years)") + 
    scale_alpha_manual(values = c(aml17 = 0.8, mds17 = 0.1))+
    theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
          legend.position = "none", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
          axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
    guides(colour = guide_legend(override.aes = list(size = rel(2))))

p <- p1+ p2 +  theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
                    legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
                    axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
                    legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
    guides(colour = guide_legend(override.aes = list(size = rel(2))), alpha="none")

leg <- cowplot::get_legend(p)
p <- p +  theme(legend.position = "none")

p <- gridExtra::grid.arrange(p1, p2, leg, ncol = 1)


dir.create(paste0(prefixOut, "/Integrated/"), showWarnings = FALSE)
ggsave(p, filename = paste0(prefixOut, "/Integrated/StateComparison.png"),
       dpi = 1600, limitsize = FALSE, height = 7, width = 6
)
