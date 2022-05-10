#/usr/bin/Rscript
library(dplyr)
library(data.table)
library(parallel)
library(dplyr)
library(Seurat)
library(Matrix)
library(future)
#pass each unbalanced training matrix files from script 1 as arg
file_name <- commandArgs(trailingOnly = T)
options(future.globals.maxSize = 1000 * 1024^2)
metadata_name <- paste0("unbalanced_training_metadata_",strsplit(file_name, "unbalanced_training_count_mat_")[[1]][2])
print(file_name)
path_to_matrix <- paste0("unbalanced_set_full_classifier/", file_name)
path_to_metadata <- paste0("unbalanced_set_full_classifier/", metadata_name)
#load training matrices (unbalanced)
count_matrix <- readRDS(path_to_matrix)
metadata <- readRDS(path_to_metadata)
colnames(count_matrix) <- row.names(metadata)
#remove ribosomal genes
ribo_genes <- grep("^RPS|^RPL", row.names(count_matrix))
count_matrix <- count_matrix[-ribo_genes,]
#randomly sample 50 000 cells since taking too long
#seurat DGE analysis
seurat_object <- CreateSeuratObject(count_matrix, meta.data = metadata)
set.seed(123)
seurat_object <- SCTransform(seurat_object)
seurat_object <- RunPCA(seurat_object, assay = "SCT")
seurat_object <- RunUMAP(seurat_object, dims = c(1:50))
seurat_object <- FindNeighbors(seurat_object, dims = 1:50)
seurat_object <- FindClusters(seurat_object, resolution = 2)

dataset_cluster_df <- seurat_object@meta.data[,c("seurat_clusters", "cell_type_merged_modified")]

#get clusters for each dataset
cell_types <- dataset_cluster_df[,"cell_type_merged_modified"]

dataset_cluster_df$cell_type_merged_modified <- cell_types

dataset_cluster <- c()
for (i in 0:(length(levels(as.factor(dataset_cluster_df$seurat_clusters)))-1)){
  dataset_cluster <- c(dataset_cluster, levels(forcats::fct_infreq(dataset_cluster_df[which(dataset_cluster_df$seurat_clusters == i),"cell_type_merged_modified"]))[1])
}
cluster_annotation <- data.table("Cluster" = 0:(length(levels(as.factor(dataset_cluster_df$seurat_clusters)))-1),
                                 "dataset" = dataset_cluster)

per_dataset_cluster_list <- list()
for (i in 1:length(levels(as.factor(cluster_annotation$dataset)))){
  per_dataset_cluster_list[[i]] <- cluster_annotation[which(cluster_annotation$dataset == levels(as.factor(cluster_annotation$dataset))[i]),1]
  names(per_dataset_cluster_list)[i] <- levels(as.factor(cluster_annotation$dataset))[i]
}
#in this script we used 20 cores for parallel processing, set workers to number of cores you would like to use
#for large layers (layer 1, layer 2 blood, layer 2 non blood) we split this for loop over multiple jobs to speed it up where we passed i as an argument
per_dataset_marker_list <- list()
for(i in 1:length(per_dataset_cluster_list)){
  all_clusters <- unlist(per_dataset_cluster_list, use.names = F)
  ident_1 <- unlist(per_dataset_cluster_list[[i]], use.names=F)
  index_remove <- c()
  for(j in 1:length(ident_1)){
    index_remove <- c(index_remove, which(all_clusters == ident_1[j]))
  }
  ident_2 <- all_clusters[-index_remove]
  plan("multiprocess", workers = 20)
  per_dataset_marker_list[[i]] <- FindMarkers(seurat_object, ident.1 = ident_1, iden.2 = ident_2, verbose = FALSE)
  print(paste0(i, " out of ",length(per_dataset_cluster_list)))
}

names(per_dataset_marker_list) <- names(per_dataset_cluster_list)
for(i in 1:length(per_dataset_cluster_list)){
  gene <- row.names(per_dataset_marker_list[[i]])
  per_dataset_marker_list[[i]] <- cbind(per_dataset_marker_list[[i]], gene)
}
#save markers list
saveRDS(per_dataset_marker_list, paste0("markers_full_classifier/markers_", strsplit(file_name, "unbalanced_training_count_mat_")[[1]][2]))
