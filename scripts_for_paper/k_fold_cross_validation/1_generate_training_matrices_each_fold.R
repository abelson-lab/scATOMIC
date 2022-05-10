#/usr/bin/Rscript

#set up matrices for layer_1
#for feature selection we use the unbalanced matrices (not class balanced)
#for training we use class balanced matrices (after k-fold split)
#layer 1
library(dplyr)
library(data.table)
library(parallel)
library(dplyr)
library(Matrix)
library(Seurat)
count_matrix <- readRDS("combined_RNA_count_matrix_layer_1.RDS")
combined_metadata <- readRDS("combined_RNA_count_matrix_layer_1_metadata.RDS")
cell_types <- combined_metadata[,"cell_type_merged"]
#remove cells that we are not including in the classifier
index_remove <- which(cell_types %in%c("Doublet", "dnT", "ILC","CD4 CTL", "gdT","MAIT","Proliferating" , "Head and Neck Cancer", "Thyroid Cancer",
                                       "Eryth", "Platelet"))
count_matrix <- count_matrix[,-index_remove]
cell_types <- cell_types[-index_remove]
#summarize cell types with higher resolution to desired resolution
cell_types <- gsub("CD14 Mono|CD16 Mono", "Macrophage", cell_types)
cell_types <- gsub("cDC1|cDC2", "cDC", cell_types)

cell_types <- gsub("B intermediate|B memory|B naive", "B cell", cell_types)
cell_types <- gsub("CD4 Naive|CD4 Proliferating|CD4 TCM|CD4 TEM|Treg", "CD4+ T cell", cell_types)
cell_types <- gsub("CD8 Naive|CD8 Proliferating|CD8 TCM|CD8 TEM", "CD8+ T cell", cell_types)
cell_types <- gsub("NK|NK Proliferating|NK_CD56bright", "Natural killer cell", cell_types)


combined_metadata <- combined_metadata[-index_remove,]

combined_metadata$cell_type_merged_modified <- cell_types

#filter dataset to include cells of interest

balanced_number <- 2500

#filter dataset
count_matrix <- count_matrix[,row.names(combined_metadata)]




#randomize dataset
set.seed(123)
randomized_index <- sample(ncol(count_matrix))
randomized_dataset <- count_matrix[,randomized_index]
randomized_metadata <- combined_metadata[randomized_index,]


set.seed(456)
#Create 5 equally size folds for each cell type!
cell_classes <- levels(as.factor(randomized_metadata$cell_type_merged_modified))
fold_per_cell_class <- list()
for(x in 1:length(cell_classes)){
  col_names_per_fold <- row.names(randomized_metadata)[which(randomized_metadata$cell_type_merged_modified == cell_classes[x])]
  fold_per_cell_class[[x]] <- cut(which(randomized_metadata$cell_type_merged_modified == cell_classes[x]),breaks=5,labels=FALSE)
  names(fold_per_cell_class[[x]]) <- col_names_per_fold

}

folds <- unlist(fold_per_cell_class)

#save matrices for each
for(j in 1:5){
  testIndexes <- which(folds==j,arr.ind=TRUE)
  testing <- randomized_dataset[, testIndexes]
  training <- randomized_dataset[,-testIndexes]
  training_metadata <- randomized_metadata[-testIndexes,]
  testing_metadata <- randomized_metadata[testIndexes,]
  colnames(training) <- row.names(training_metadata)
  training_seurat <- CreateSeuratObject(training)
  nFeatureRNA <- training_seurat@meta.data$nFeature_RNA
  names(nFeatureRNA) <- row.names(training_seurat@meta.data)
  rm(training_seurat)
  training <- training[, intersect(names(which(nFeatureRNA > 500 & nFeatureRNA < 6000)), colnames(training))]
  training_metadata <- training_metadata[colnames(training), ]
  saveRDS(training, paste0("/unbalanced_set/unbalanced_training_count_mat_layer_1_fold_", j, ".RDS"))
  saveRDS(training_metadata, paste0("/unbalanced_set/unbalanced_training_metadata_layer_1_fold_", j, ".RDS"))
  saveRDS(testing, paste0("/testing_sets/unbalanced_testing_count_mat_layer_1_fold_", j, ".RDS"))
  saveRDS(testing_metadata, paste0("/testing_sets/unbalanced_testing_metadata_layer_1_fold_", j, ".RDS"))

  #sample each class
  cell_types_merged <- levels(as.factor(as.character(training_metadata$cell_type_merged_modified)))
  index_keep <- c()
  for(i in 1:length(cell_types_merged)){
    index_cell_type <- which(training_metadata$cell_type_merged_modified == cell_types_merged[i])
    if(length(index_cell_type) >= balanced_number){
      set.seed(456)
      index_keep <- c(index_keep, sample(index_cell_type, balanced_number))
    } else if (length(index_cell_type) < balanced_number){
      index_keep <- c(index_keep, sample(index_cell_type, balanced_number, replace = T))
    }
  }
  training_metadata$cell_type_merged_modified <- as.factor(as.character(training_metadata$cell_type_merged_modified))
  training_metadata <- training_metadata[index_keep,]
  training <- training[,index_keep]

  saveRDS(training, paste0("/split_by_folds/balanced_training_count_mat_layer_1_fold_", j, ".RDS"))
  saveRDS(training_metadata, paste0("/split_by_folds/balanced_training_metadata_layer_1_fold_", j, ".RDS"))
  print(j)
}

