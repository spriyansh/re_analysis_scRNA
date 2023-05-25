# Call the required libraries
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))

# Prefix
prefixIn <- "../data/input/  "
prefixOut <- "../data/output/"
prefixMDS_Genes <- "../data/output/MDS_genes/"
prefixTrMaps <- "../data/output/tr_maps/"
prefixTrend <- "../data/output/trend/"

# mds genes
mds_genes <- c(
    "TET2", "SF3B1", "ASXL1", "DNMT3A", "SRSF2", "RUNX1", "TP53", "U2AF1",
    "EZH2", "ZRSR2", "STAG2", "CBL", "NRAS", "JAK2", "SETBP1", "IDH1", "IDH2", "ETV6", "FLT3",
    "NF1", "CALR", "MPL", "GATA2"
)

# Other Explorative Genes
oth_genes <- c("MPO", "CD34", "CD79B", "GATA1", "IRF8", "CD41",  "EPOR", "EPOR","AICDA","APOBEC3C", "APOBEC3F","APOBEC3G")

# all_file_names
all_file_names <- list.files(paste0(prefixOut))

# Distribute filenames
azi_obj_name <- grep(all_file_names, pattern = "_final_obj.RDS", value = T)

trash <- capture.output(
# Run for every dataset
for (i in azi_obj_name) {
    
    # Get the Replicate
    #i <- azi_obj_name[1]
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
        age = 35
        gen = "Male"
        
    }else if (names(rep_split) == "rep2"){
        pca_dim = 5
        umap_min_dist = 0.9
        umap.n_neighbors = 25
        age = 28
        gen = "Female"
    }else if (names(rep_split) == "rep3"){
        pca_dim = 5
        umap_min_dist = 0.8
        umap.n_neighbors = 15
        age = 19
        gen = "Female"
    }
    
    # Monocel3 Processing
    cds <- preprocess_cds(cds, num_dim = pca_dim, method = "PCA", norm_method = "log")
    cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method="UMAP",
                            umap.min_dist = umap_min_dist, umap.n_neighbors= umap.n_neighbors)
    cds <- cluster_cells(cds, reduction_method = "UMAP")
    
    # Plotting
    png(paste0(prefixTrMaps, names(rep_split), "_UMAP_Cells.png"), width = 1000, height = 1000)
    p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 2, 
               color_cells_by = "az_labels", label_cell_groups = T)+ 
        scale_color_brewer(palette = "Dark2") + 
        ggtitle(paste("Replicate:",names(rep_split), "|","Age:", age, "|", " Gender:", gen)) +
    theme(legend.position = "bottom", axis.text = element_text(size = rel(3)), legend.title = element_blank(),
            axis.title = element_text(size = rel(4)), legend.text = element_text(size = rel(3)),
          legend.key.size = unit(rel(2), 'cm'), plot.title = element_text(size = rel(4))) + 
        guides(colour = guide_legend(override.aes = list(size=10)))
    print(p)
    dev.off()
    
    # Learn Graph
    cds <- learn_graph(cds)
    
    # Plotting
    png(paste0(prefixTrMaps, names(rep_split), "_UMAP_TLine.png"), width = 1000, height = 1000)
    p1 <- plot_cells(cds, reduction_method = "UMAP", cell_size = 2, 
               color_cells_by = "az_labels", label_cell_groups = T, trajectory_graph_segment_size = 3)+ 
        scale_color_brewer(palette = "Dark2") + 
        ggtitle(paste("Replicate:",names(rep_split), "|","Age:", age, "|", " Gender:", gen))+
        theme(legend.position = "bottom", axis.text = element_text(size = rel(3)),
              legend.title = element_blank(), axis.title = element_text(size = rel(4)),
              legend.text = element_text(size = rel(3)), legend.key.size = unit(rel(2), 'cm'),
              plot.title = element_text(size = rel(4))) +
        guides(colour = guide_legend(override.aes = list(size=10)))
    print(p1)
    dev.off()
    
    # Order cells
    cds <- order_cells(cds, root_cells = rownames(colData(cds)[colData(cds)$az_labels == "HSC",]))
    
    # Plotting
    png(paste0(prefixTrMaps, names(rep_split), "_UMAP_Trajectory.png"), width = 1000, height = 1000)
    p2 <- plot_cells(cds, reduction_method = "UMAP", cell_size = 2,  trajectory_graph_segment_size = 3,
               color_cells_by = "pseudotime", label_cell_groups = T) + 
        scale_color_viridis(option = "E", direction = -1)+
        ggtitle(paste("Replicate:",names(rep_split), "|","Age:", age, "|", " Gender:", gen))+ 
         theme(legend.position = "bottom", axis.text = element_text(size = rel(3)),
               legend.title = element_blank(), axis.title = element_text(size = rel(4)),
               legend.text = element_text(size = rel(3)), legend.key.size = unit(rel(2), 'cm'),
               plot.title = element_text(size = rel(4))) #+ guides(colour = guide_legend(override.aes = list(size=10)))
    print(p2)
    dev.off()
    
    png(paste0(prefixTrMaps, names(rep_split), "_UMAP_Combined.png"), width = 2000, height = 1250)
    p3 <- ggpubr::ggarrange(p1,p2)
    print(p3)
    dev.off()
    
    # Plot MDS Genes
    for(j in mds_genes){
        png(paste0(prefixMDS_Genes, names(rep_split),"_", j, "_UMAP.png"), width = 1200, height = 1200)
        p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 2, 
                        label_cell_groups = F, show_trajectory_graph = F,
                        genes = j) + 
            scale_color_viridis(option = "A", direction = -1)+
            ggtitle(paste("Gene:" ,j,"Replicate:",names(rep_split), "|","Age:", age, "|", " Gender:", gen))+ 
            theme(legend.position = "bottom", 
                  axis.text = element_text(size = rel(3)),
                  legend.title = element_blank(), axis.title = element_text(size = rel(4)),
                  legend.text = element_text(size = rel(3)), legend.key.size = unit(rel(2), 'cm'),
                  plot.title = element_text(size = rel(4)))
        print(p)
        dev.off()
    }
    
    oth_genes <- rownames(cds)[rownames(cds) %in% oth_genes]
    
    # Plot MDS Genes
    for(j in oth_genes){
        png(paste0(prefixMDS_Genes, names(rep_split),"_", j, "_UMAP.png"), width = 1200, height = 1200)
        p <- plot_cells(cds, reduction_method = "UMAP", cell_size = 2, 
                        label_cell_groups = F, show_trajectory_graph = F,
                        genes = j) + 
            scale_color_viridis(option = "A", direction = -1)+
            ggtitle(paste("Gene:" ,j,"Replicate:",names(rep_split), "|","Age:", age, "|", " Gender:", gen))+ 
            theme(legend.position = "bottom", 
                  axis.text = element_text(size = rel(3)),
                  legend.title = element_blank(), axis.title = element_text(size = rel(4)),
                  legend.text = element_text(size = rel(3)), legend.key.size = unit(rel(2), 'cm'),
                  plot.title = element_text(size = rel(4)))
        print(p)
        dev.off()
    }
    
    
    for(j in c(mds_genes,oth_genes)){ 
    
    png(paste0(prefixTrend, names(rep_split), "_",j, "_Marker_XY.png"), height = 700, width = 1300)
        cds_sub <- cds[rowData(cds)$gene_short_name %in% j,]
    p <- plot_genes_in_pseudotime(cds_sub, color_cells_by = "az_labels",cell_size = 3,
                                  label_by_short_name = F,
                             min_expr = 0.05) +  scale_color_brewer(palette = "Dark2")+
        ggtitle(paste("Gene:" ,j,"Replicate:",names(rep_split), "|","Age:", age, "|", " Gender:", gen))+
        theme(legend.position = "bottom", 
              axis.text = element_text(size = rel(3)),
              legend.title = element_blank(), axis.title = element_text(size = rel(4)),
              legend.text = element_text(size = rel(3)), legend.key.size = unit(rel(2), 'cm'),
              plot.title = element_text(size = rel(4))) +
        guides(colour = guide_legend(override.aes = list(size=10)))
     
    print(p)
    dev.off()
    }
    
    # Save Data 
    rds_name <- paste0(prefixOut, names(rep_split), "_cds_obj.RDS")
    saveRDS(cds, rds_name)
    
    cat(paste0("\nCompleted for ", rep_split[[1]][1], "\n"))
}
)
