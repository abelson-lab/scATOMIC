#' get_new_scATOMIC_layer
#'
#' @param training_data count matrix, preferably sparse matrix with rows being genes and columns being cells
#' @param cell_type_metadata vector containing the cell type corresponding to each cell in training_data
#' @param output_dir directory to save random forest and selected genes to
#' @param layer_name name of resulting layer
#' @param n_cells_replicate number of cells for balancing each class, we suggest 10000 if <5 cell types, 5000 if < 20 cell types and 2500 if > 20
#' @param n_trees number of trees in random forest
#' @param custom_gene_list vector of features to use, if user wants to provide their own genes.
#'
#' @return list of new rf model and selected features
#' @export
#' @examples
#' \dontrun{
#' #visit github for more information on how to generate training datasets.
#' breast_cancer_subclassification <- get_new_scATOMIC_layer(training_data = Wu_et_al_breast_count_mat,cell_type_metadata = Wu_et_al_breast_metadata$subtype,
#' output_dir = "Wu_etal_2021_BRCA_scRNASeq/breast_subclassification_layer/",
#'
#'
#'
#'
#' }

get_new_scATOMIC_layer <- function(training_data, cell_type_metadata,output_dir,layer_name,n_cells_replicate = 5000, n_trees = 500, custom_gene_list = NULL){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  print("Training new layers can take a long time...")
  library(dplyr)
  cell_type_metadata <- data.frame(cell_type_metadata)
  colnames(cell_type_metadata) <- "cell_type"
  row.names(cell_type_metadata) <- colnames(training_data)
  if(length(custom_gene_list) == 0){
    print("Creating Seurat Object")
    seurat_object <- Seurat::CreateSeuratObject(training_data, meta.data = cell_type_metadata)
    seurat_object <- Seurat::SCTransform(seurat_object)
    seurat_object <- Seurat::RunPCA(seurat_object, assay = "SCT")
    seurat_object <- Seurat::RunUMAP(seurat_object, dims = c(1:50))
    seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:50)
    seurat_object <- Seurat::FindClusters(seurat_object, resolution = 2)

    dataset_cluster_df <- seurat_object@meta.data[,c("seurat_clusters", "cell_type")]

    #get clusters for each dataset



    cell_types <- dataset_cluster_df[,"cell_type"]



    dataset_cluster_df$cell_type <- cell_types





    dataset_cluster <- c()
    for (i in 0:(length(levels(as.factor(dataset_cluster_df$seurat_clusters)))-1)){
      dataset_cluster <- c(dataset_cluster, levels(forcats::fct_infreq(dataset_cluster_df[which(dataset_cluster_df$seurat_clusters == i),"cell_type"]))[1])
    }
    cluster_annotation <- data.table::data.table("Cluster" = 0:(length(levels(as.factor(dataset_cluster_df$seurat_clusters)))-1),
                                                 "dataset" = dataset_cluster)

    per_dataset_cluster_list <- list()
    for (i in 1:length(levels(as.factor(cluster_annotation$dataset)))){
      per_dataset_cluster_list[[i]] <- cluster_annotation[which(cluster_annotation$dataset == levels(as.factor(cluster_annotation$dataset))[i]),1]
      names(per_dataset_cluster_list)[i] <- levels(as.factor(cluster_annotation$dataset))[i]
    }
    print("running DGE analysis")
    per_dataset_marker_list <- list()
    for(i in 1:length(per_dataset_cluster_list)){
      all_clusters <- unlist(per_dataset_cluster_list, use.names = F)
      ident_1 <- unlist(per_dataset_cluster_list[[i]], use.names=F)
      index_remove <- c()
      for(j in 1:length(ident_1)){
        index_remove <- c(index_remove, which(all_clusters == ident_1[j]))
      }
      ident_2 <- all_clusters[-index_remove]
      per_dataset_marker_list[[i]] <- Seurat::FindMarkers(seurat_object, ident.1 = ident_1, iden.2 = ident_2)

    }
    names(per_dataset_marker_list) <- names(per_dataset_cluster_list)
    for(i in 1:length(per_dataset_cluster_list)){
      gene <- row.names(per_dataset_marker_list[[i]])
      per_dataset_marker_list[[i]] <- cbind(per_dataset_marker_list[[i]], gene)
    }
    print("Training classifier")
    number_sd <- 0
    percent_cutoff <- 0.4
    cell_types <- levels(as.factor(as.character(cell_type_metadata$cell_type)))
    index_keep <- c()
    for(i in 1:length(cell_types)){
      index_cell_type <- which(cell_type_metadata$cell_type == cell_types[i])
      if(length(index_cell_type) >= n_cells_replicate){
        set.seed(456)
        index_keep <- c(index_keep, sample(index_cell_type, n_cells_replicate))
      } else if (length(index_cell_type) < n_cells_replicate){
        index_keep <- c(index_keep, sample(index_cell_type, n_cells_replicate, replace = T))
      }
    }
    cell_type_metadata$cell_type <- as.factor(as.character(cell_type_metadata$cell_type))
    cell_type_metadata <- cell_type_metadata[index_keep,]
    training_data <- training_data[,index_keep]
    colnames(training_data) <- row.names(cell_type_metadata)


    for (j in 1:length(per_dataset_marker_list)){
      diff <- per_dataset_marker_list[[j]]$pct.1 - per_dataset_marker_list[[j]]$pct.2
      per_dataset_marker_list[[j]] <- cbind(per_dataset_marker_list[[j]], diff)
      diffxFC <- per_dataset_marker_list[[j]]$avg_log2FC * per_dataset_marker_list[[j]]$diff
      per_dataset_marker_list[[j]] <- cbind(per_dataset_marker_list[[j]], diffxFC)
      per_dataset_marker_list[[j]] <- per_dataset_marker_list[[j]][which(per_dataset_marker_list[[j]]$diff > 0),]
      per_dataset_marker_list[[j]] <- per_dataset_marker_list[[j]][which(per_dataset_marker_list[[j]]$pct.2 < percent_cutoff),]

    }

    top_pos_per_test <- list()
    #set n to the number of features per cluster to extract
    for (m in 1:length(per_dataset_marker_list)){
      cutoff_per_class = mean(per_dataset_marker_list[[m]]$diff) + number_sd*sd(per_dataset_marker_list[[m]]$diff)
      top_pos_per_test[[m]] <- per_dataset_marker_list[[m]][which(per_dataset_marker_list[[m]]$diff > cutoff_per_class),]
      if(nrow(top_pos_per_test[[m]]) < 50){
        top_pos_per_test[[m]] <- per_dataset_marker_list[[m]] %>% dplyr::top_n(n = 50, wt = diff)
      }
      if(nrow(top_pos_per_test[[m]]) > 200){
        top_pos_per_test[[m]] <- top_pos_per_test[[m]] %>% dplyr::top_n(n = 200, wt = diff)
      }


      top_pos_per_test[[m]] <- top_pos_per_test[[m]]$gene
    }
    names(top_pos_per_test)  <- names(per_dataset_marker_list)
    top_pos_per_test <- as.vector(unlist(top_pos_per_test, use.names = FALSE))
    top_genes_unlisted <- unique(top_pos_per_test)

  }else{
    print("Training classifier")
    cell_types <- levels(as.factor(as.character(cell_type_metadata$cell_type)))
    index_keep <- c()
    for(i in 1:length(cell_types)){
      index_cell_type <- which(cell_type_metadata$cell_type == cell_types[i])
      if(length(index_cell_type) >= n_cells_replicate){
        set.seed(456)
        index_keep <- c(index_keep, sample(index_cell_type, n_cells_replicate))
      } else if (length(index_cell_type) < n_cells_replicate){
        index_keep <- c(index_keep, sample(index_cell_type, n_cells_replicate, replace = T))
      }
    }
    cell_type_metadata$cell_type <- as.factor(as.character(cell_type_metadata$cell_type))
    cell_type_metadata <- cell_type_metadata[index_keep,]
    training_data <- training_data[,index_keep]
    colnames(training_data) <- row.names(cell_type_metadata)


    top_genes_unlisted <- unique(custom_gene_list)
  }

  training_data <- as(t(as.matrix(training_data)), "sparseMatrix")
  training_cell_lines_pbmc_normalized_counts <- Rmagic::library.size.normalize(training_data)
  training_cell_lines_pbmc_normalized_counts <-as(t(as.matrix(training_cell_lines_pbmc_normalized_counts)), "sparseMatrix")
  training_cell_lines_pbmc_norm_transformed_counts_filtered_combo <- na.omit(training_cell_lines_pbmc_normalized_counts[top_genes_unlisted,])
  training_scrna_counts_as_fraction <- as.data.frame(apply(training_cell_lines_pbmc_norm_transformed_counts_filtered_combo, 2, function(x){
    x/sum(x)
  }))
  training_scaled_fractions_scRNA <- as.data.frame(t(training_scrna_counts_as_fraction))
  cell_class <- cell_type_metadata
  training_scaled_fractions_scRNA <- cbind(cell_class, training_scaled_fractions_scRNA)
  names_genes <- colnames(training_scaled_fractions_scRNA)
  names_genes <- gsub("-", "_", names_genes)
  colnames(training_scaled_fractions_scRNA) <- names_genes
  training_scaled_fractions_scRNA <- na.omit(training_scaled_fractions_scRNA)
  print("training started")

  library(randomForest)

  rf_classifier_cell_lines <- randomForest::randomForest(cell_class ~ ., data=training_scaled_fractions_scRNA, ntree = n_trees)



  print("training ended")
  # keep genes expressed in at least 10 cells
  print("saving")
  rm(list = setdiff(ls(), c("rf_classifier_cell_lines", "output_dir", "layer_name",
                            "top_genes_unlisted")))
  stripRF <- function(cm) {
    cm$finalModel$predicted <- NULL
    cm$finalModel$oob.times <- NULL
    cm$finalModel$y <- NULL
    cm$finalModel$votes <- NULL
    cm$control$indexOut <- NULL
    cm$control$index    <- NULL
    cm$trainingData <- NULL

    attr(cm$terms,".Environment") <- c()
    attr(cm$formula,".Environment") <- c()

    cm
  }
  rf_classifier_cell_lines <- stripRF(rf_classifier_cell_lines)
  saveRDS(rf_classifier_cell_lines, paste0(output_dir,layer_name, "_random_forest_classifier.RDS"))
  saveRDS(top_genes_unlisted, paste0(output_dir,layer_name, "_selected_features.RDS"))

  print(paste0(layer_name, " saved to ", output_dir))
  new_model <- list()
  new_model[[1]] <- rf_classifier_cell_lines
  new_model[[2]] <- top_genes_unlisted
  names(new_model) <- c("rf_classifier", "selected_features")
  return(new_model)
}
