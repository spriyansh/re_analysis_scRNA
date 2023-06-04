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

# Prefix
prefixIn <- "guess_et_al/data/output/p17/Azimuth/"
dir.create("guess_et_al/data/output/p17/Monocel3/", showWarnings = FALSE)
prefixOut <- "guess_et_al/data/output/p17/Monocel3/"

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
all_file_names <- grep(all_file_names, pattern = ".h5seurat", value = T)

# Distribute filenames
aml17_names <- grep(all_file_names, pattern = "aml17", value = T) # Male-35
mds17_names <- grep(all_file_names, pattern = "mds17", value = T) # Female-28

# Make a list of filenames
rep_list <- list(aml17 = aml17_names, 
                 mds17 = mds17_names)

  # Run for every dataset
  for (i in names(rep_list)) {
      
      print(i)
      if (i == "mds17") {
          patientID <- "17"
          condition <- "MDS"
          umap_min_dist <- 0.1
          umap.n_neighbors <- 15
          wd =15
          res = 0.002
      } else if (i == "aml17") {
          patientID <- "17"
          condition <- "sAML"
          umap_min_dist <- 0.1
          umap.n_neighbors <- 15
          wd =15
          res = 0.002
  } 
      
      # Get the names rep#
      rep_i_name <- rep_list[[i]]
      
      #  the Processed Seurat Object
      sob <- LoadH5Seurat(file = paste0(prefixIn, rep_i_name), verbose = F)
      
      # Create Monocel3 Object
      cds <- new_cell_data_set(
          expression_data = sob@assays$RNA@data,
          cell_metadata = sob@meta.data,
          gene_metadata = data.frame(
              gene_short_name = rownames(sob@assays$RNA@data),
              row.names = rownames(sob@assays$RNA@data)
          )
      )

    # Monocel3 Processing
    cds <- preprocess_cds(cds, num_dim = 5, method = "PCA", norm_method = "none")
    cds <- reduce_dimension(cds,
      preprocess_method = "PCA", reduction_method = "UMAP",
      umap.min_dist = umap_min_dist, umap.n_neighbors = umap.n_neighbors
    )
    
    # Plotting Monocle3 UMAP
    p <- plot_cells(cds,
               reduction_method = "UMAP", cell_size = 0.5,
               color_cells_by = "cell.type", label_cell_groups = F
    ) +
        #scale_color_brewer(palette = "Set1") +
        ggtitle(paste0("PatientID: ", patientID," Condition: ", condition)) +
        theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
            legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
            axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
            legend.key.size = unit(0.1, "cm"), plot.title = element_text(size = rel(1))) +
        guides(colour = guide_legend(override.aes = list(size = rel(4))))
    dir.create(paste0(prefixOut, "/UMAP"), showWarnings = FALSE)
    ggsave(p, filename = paste0(prefixOut, "/UMAP/", i, "_UMAP.png"),
           dpi = 1300, limitsize = FALSE,width = 7
    )
    
    # Perform Liden Clustering
    cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = res)
    
    # Overring Partition
    # partition_vec <- rep(1, ncol(cds))
    # names(partition_vec) <- colnames(cds)
    # cds@clusters@listData[["UMAP"]][["partitions"]] <- as.factor(partition_vec)

    # Plotting Monocle3 UMAP
    p <- plot_cells(cds,
                    reduction_method = "UMAP", cell_size = 0.5,
                    color_cells_by = "cluster", label_cell_groups = F
    ) +
        scale_color_brewer(palette = "Set1") +
        ggtitle(paste0("PatientID: ", patientID," Condition: ", condition)) +
        theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
              legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
              axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
              legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
        guides(colour = guide_legend(override.aes = list(size = rel(6))))
    dir.create(paste0(prefixOut, "/Cluster"), showWarnings = FALSE)
    ggsave(p, filename = paste0(prefixOut, "/Cluster/", i, "_Cluster.png"),
           dpi = 1300, limitsize = FALSE
    )

    # Learn Graph
    cds <- learn_graph(cds,verbose = F)
    
    # Plotting Monocle3 UMAP
    p1 <- plot_cells(cds,
                    reduction_method = "UMAP", cell_size = 0.5,
                    color_cells_by = "cell.type", label_cell_groups = F
    ) +
        #scale_color_brewer(palette = "Dark2") +
        ggtitle(paste0("PatientID: ", patientID," Condition: ", condition)) +
        theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
              legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
              axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
              legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) +
        guides(colour = guide_legend(override.aes = list(size = rel(6))))
    dir.create(paste0(prefixOut, "/TrajectoryLine"), showWarnings = FALSE)
    ggsave(p1, filename = paste0(prefixOut, "/TrajectoryLine/", i, "_TrLine.png"),
           dpi = 1300, limitsize = FALSE
    )
    
    # Order cells
    cds <- order_cells(cds, root_cells = rownames(
        colData(cds)[colData(cds)$cell.type == "HSC", ]))

    # Plotting Monocle3 UMAP
    p2 <- plot_cells(cds,
                    reduction_method = "UMAP", cell_size = 0.5,
                    color_cells_by = "pseudotime", label_cell_groups = F,
                    label_branch_points = F,label_leaves = F,
                    label_roots = T, label_principal_points = F
    ) +
        scale_color_viridis(option = "C") +
        ggtitle(paste0("PatientID: ", patientID," Condition: ", condition)) +
        theme(plot.margin = margin(0.2,0.2,0,0, "cm"),
              legend.position = "bottom", axis.text = element_text(size = rel(1)), legend.title = element_blank(),
              axis.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)),
              legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = rel(1))) 
    dir.create(paste0(prefixOut, "/Pseudotime"), showWarnings = FALSE)
    ggsave(p2, filename = paste0(prefixOut, "/Pseudotime/", i, "_Pseudotime.png"),
           dpi = 1300, limitsize = FALSE
    )
   
    p <- p1+p2 
    ggsave(p, filename = paste0(prefixOut, "/Pseudotime/", i, "_PseudotimeCombined.png"),
           dpi = 1300, limitsize = FALSE, width = wd
    )
    
    # Save the Object
    saveRDS(cds, paste0(prefixOut, i, "_m3_UMAP_OB.RDS"))
    
    cat(paste0("\nCompleted for ", i, "\n"))
  }
