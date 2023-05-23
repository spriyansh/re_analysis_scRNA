# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(RColorBrewer))

# Prefix
prefixIn <- "../input/"
prefixOut <- "../output/"
prefixMDS_Genes <- "../output/MDS_genes/"
prefixTrMaps <- "../output/tr_maps/"

# mds genes
mds_genes <- c(
    "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
    "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
    "NF1", "CALR", "MPL", "GATA2"
)

# Other Explorative Genes
oth_genes <- c("MPO", "CD34", "CD79B", "GATA1", "IRF8", "CD41")

# all_file_names
all_file_names <- list.files(paste0(prefixOut))

# Distribute filenames
azi_obj_name <- grep(all_file_names, pattern = "_final_obj.RDS", value = T)

# Run for every dataset
for (i in azi_obj_name) {
    
    # Get the Replicate
    rep_split <- str_split(i, "_")
    names(rep_split) <- rep_split[[1]][1]
    
    # Read the dataset
    s.obj <- readRDS(paste0(prefixOut, i))
    
    # Create cds
    cds <- new_cell_data_set(expression_data = s.obj@assays$RNA@counts,
                             cell_metadata = s.obj@meta.data,
                             gene_metadata = data.frame(
                                 gene_short_name = rownames(s.obj@assays$RNA@counts),
                                 row.names = rownames(s.obj@assays$RNA@counts)
                             )
    )
    
    if (names(rep_split) == "rep1"){
        pca_dim = 5
        umap_min_dist = 0.8
        umap.n_neighbors = 15
    }else if (names(rep_split) == "rep2"){
        pca_dim = 5
        umap_min_dist = 0.9
        umap.n_neighbors = 25
    }else if (names(rep_split) == "rep2"){
        pca_dim = 5
        umap_min_dist = 0.8
        umap.n_neighbors = 15
    }
    
    # Monocel3 Processing
    cds <- preprocess_cds(cds, num_dim = pca_dim, method = "PCA", norm_method = "log")
    cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method="UMAP",
                            umap.min_dist = umap_min_dist, umap.n_neighbors= umap.n_neighbors)
    cds <- cluster_cells(cds, reduction_method = "UMAP")
    
    # Plotting
    png(paste0(prefixTrMaps, names(rep_split), "UMAP_Cells.png"))
    p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 1.3, 
               color_cells_by = "az_labels", label_cell_groups = T)+ 
        scale_color_brewer(palette = "Dark2") + theme(legend.position = "bottom")
    print(p)
    dev.off()
    
    # Learn Graph
    cds <- learn_graph(cds)
    
    # Plotting
    png(paste0(prefixTrMaps, names(rep_split), "UMAP_TLine.png"))
    p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 1.3, 
               color_cells_by = "az_labels", label_cell_groups = T)+ 
        scale_color_brewer(palette = "Dark2") + theme(legend.position = "bottom")
    print(p)
    dev.off()
    
    # Order cells
    cds <- order_cells(cds, root_cells = rownames(colData(cds)[colData(cds)$az_labels == "HSC",]))
    
    # Plotting
    png(paste0(prefixTrMaps, names(rep_split), "UMAP_Trajectory.png"))
    p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 1.3, 
               color_cells_by = "pseudotime", label_cell_groups = T)+ 
         theme(legend.position = "bottom")
    print(p)
    dev.off()
    
    # Plot MDS Genes
    for(j in mds_genes){
        png(paste0(prefixMDS_Genes, names(rep_split),"_", j, "_UMAP.png"))
        p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 1.3, 
                        label_cell_groups = F, show_trajectory_graph = F,
                        genes = j)+ 
            theme(legend.position = "bottom")
        print(p)
        dev.off()
    }
    
    oth_genes <- rownames(cds)[rownames(cds) %in% oth_genes]
    
    # Plot MDS Genes
    for(j in oth_genes){
        png(paste0(prefixMDS_Genes, names(rep_split),"_", j, "_UMAP.png"))
        p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 1.3, 
                        label_cell_groups = F, show_trajectory_graph = F,
                        genes = c("APOBEC3C", "APOBEC3F","APOBEC3G","EPOR"))+ 
            theme(legend.position = "bottom")
        print(p)
        dev.off()
    }
    
    cds_sub <- cds[rowData(cds)$gene_short_name %in% c("EPOR","AICDA","APOBEC3C", "APOBEC3F","APOBEC3G"),]
    
    
    png(paste0(prefixMDS_Genes, names(rep_split), "Marker_XY.png"), height = 800, width = 500)
    plot_genes_in_pseudotime(cds_sub,
                             color_cells_by = "az_labels",
                             min_expr = 0.05) +  scale_color_brewer(palette = "Dark2") 
    print(p)
    dev.off()
    
    # Save Data 
    rds_name <- paste0(prefixOut, names(rep_split), "_cds_obj.RDS")
    saveRDS(cds, rds_name)
    
    cat(paste0("\nCompleted for ", rep_split[[1]][1], "\n"))
}
