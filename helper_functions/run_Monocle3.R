suppressPackageStartupMessages(require(monocle3))

run_monocel3 <- function(object, object.type = "seurat",
                         umap.dim = 2, pca.dim = 10,
                         cluster.resolution = NULL,
                         umap_min_dist = 0.3,
                         umap.n_neighbors = 20,
                         norm_method = "none",
                         progen.col = "Azimuth.Cells",
                         progen.start = "HSC"){
    
    if (object.type == "seurat"){
        # Create Monocel3 Object
        cds <- new_cell_data_set(
            expression_data = object@assays$RNA@data,
            cell_metadata = object@meta.data,
            gene_metadata = data.frame(
                gene_short_name = rownames(object@assays$RNA@data),
                row.names = rownames(object@assays$RNA@data)
                )
        )
    }else if(object.type == "monocle3"){
        cds <- object
    }else{
        print("Please use 'Seurat' for now")
    }
    
    # Preprocess
    cds <- preprocess_cds(cds, num_dim = pca.dim, method = "PCA", norm_method = norm_method)
    
    # Reduce Dim
    cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP",
                            max_components = umap.dim, umap.min_dist = umap_min_dist, 
                            umap.n_neighbors = umap.n_neighbors
                            )
    
    # Perform Liden Clustering
    cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = cluster.resolution)
    
    # Learn Graph
    cds <- learn_graph(cds, verbose = F)
    
    # Order cells
    cds <- order_cells(cds, root_cells = rownames(
        colData(cds)[colData(cds)[[progen.col]] == progen.start, ]
    ))
    
    # Return
    return(cds)
}