#/usr/bin/Rscript

#set up matrices for training each classification branch
#for feature selection we use the raw matrices (not class balanced)
#for training we use class balanced matrices (after k-fold split)
#layer 1
library(dplyr)
library(data.table)
library(parallel)
library(dplyr)
library(Matrix)
library(Seurat)
#k fold index (1:5)
index_for_loop <- as.numeric(commandArgs(trailingOnly = T))[1]
layer_name <- as.numeric(commandArgs(trailingOnly = T))[2]

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
if(layer_name == "layer_1"){
  cells_keep <- levels(as.factor(cell_types))
  balanced_number <- 2500
} else if(layer_name == "layer_2_non_blood"){
  cells_keep <- c("Bile Duct Cancer","Bladder Cancer","Bone Cancer","Brain Cancer",
                  "Breast Cancer","Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Endothelial Cells",
                  "Esophageal Cancer","Fibroblasts","Gallbladder Cancer","Gastric Cancer",
                  "Glial Cells","Kidney Cancer","Liver Cancer","Lung Cancer",
                  "Myofibroblasts","Neuroblastoma","Oligodendrocytes","Ovarian Cancer",
                  "Pancreatic Cancer","Prostate Cancer","Rhabdomyosarcoma","Sarcoma",
                  "Skin Cancer","Smooth Muscle Cells")
  balanced_number <- 5000

} else if(layer_name == "layer_2_blood"){
  cells_keep <- c("B cell","CD4+ T cell",
                  "CD8+ T cell", "HSPC",
                  "Macrophage","Natural killer cell", "Plasmablast", "cDC", "pDC", "ASDC")
  balanced_number <- 10000

} else if(layer_name == "layer_3_B_cell"){
  cells_keep <- c("B cell", "Plasmablast")
  balanced_number <- 10000

} else if(layer_name == "layer_3_MDC"){
  cells_keep <- c("cDC", "pDC", "ASDC", "Macrophage")
  balanced_number <- 10000

} else if(layer_name == "layer_3_non_stromal"){
  cells_keep <- c("Bile Duct Cancer","Bladder Cancer","Bone Cancer","Brain Cancer",
                  "Breast Cancer","Colon/Colorectal Cancer","Endometrial/Uterine Cancer",
                  "Esophageal Cancer","Gallbladder Cancer","Gastric Cancer",
                  "Kidney Cancer","Liver Cancer","Lung Cancer",
                  "Neuroblastoma","Ovarian Cancer",
                  "Pancreatic Cancer","Prostate Cancer","Rhabdomyosarcoma","Sarcoma",
                  "Skin Cancer")
  balanced_number <- 5000

} else if(layer_name == "layer_3_stromal"){
  cells_keep <- c("Endothelial Cells",
                  "Fibroblasts",
                  "Myofibroblasts","Smooth Muscle Cells")
  balanced_number <- 5000

} else if(layer_name == "layer_3_TNK"){
  cells_keep <- c("CD4+ T cell",
                  "CD8+ T cell", "Natural killer cell")
  balanced_number <- 10000

} else if(layer_name == "layer_4_CD4_CD8"){
  cells_keep <- c("CD4+ T cell",
                  "CD8+ T cell")
  balanced_number <- 10000

} else if(layer_name == "layer_4_CD8_NK"){
  cells_keep <- c("Natural killer cell",
                  "CD8+ T cell")
  balanced_number <- 10000

} else if(layer_name == "layer_4_DC"){
  cells_keep <- c("cDC","pDC",
                  "ASDC")
  balanced_number <- 10000

} else if(layer_name == "layer_4_non_GI"){
  cells_keep <- c("Breast Cancer","Endometrial/Uterine Cancer",
                  "Kidney Cancer","Lung Cancer",
                  "Ovarian Cancer",
                  "Prostate Cancer")
  balanced_number <- 5000

} else if(layer_name == "layer_4_GI"){
  cells_keep <- c("Bile Duct Cancer","Bladder Cancer","Colon/Colorectal Cancer",
                  "Esophageal Cancer","Gallbladder Cancer","Gastric Cancer",
                  "Liver Cancer","Pancreatic Cancer")
  balanced_number <- 5000

} else if(layer_name == "layer_4_soft_tissue_neuro"){
  cells_keep <- c("Bone Cancer","Brain Cancer","Neuroblastoma",
                  "Lung Cancer", "Rhabdomyosarcoma","Sarcoma",
                  "Skin Cancer")
  balanced_number <- 5000

} else if(layer_name == "layer_5_soft_tissue_neuro"){
  cells_keep <- c("Bone Cancer","Brain Cancer",
                  "Neuroblastoma","Rhabdomyosarcoma","Sarcoma")
  balanced_number <- 5000

} else if(layer_name == "layer_5_breast_lung_prostate"){
  cells_keep <- c("Breast Cancer", "Lung Cancer", "Prostate Cancer")
  balanced_number <- 10000

} else if(layer_name == "layer_5_ov_endo_kid"){
  cells_keep <- c("Endometrial/Uterine Cancer","Kidney Cancer","Ovarian Cancer")
  balanced_number <- 10000

} else if(layer_name == "layer_5_digestive"){
  cells_keep <- c("Colon/Colorectal Cancer",
                  "Esophageal Cancer","Gastric Cancer")
  balanced_number <- 10000

} else if(layer_name == "layer_5_biliary"){
  cells_keep <- c("Bile Duct Cancer","Bladder Cancer", "Gallbladder Cancer",
                  "Liver Cancer","Pancreatic Cancer")
  balanced_number <- 10000

} else if(layer_name == "layer_6_soft_tissue"){
  cells_keep <- c("Bone Cancer","Rhabdomyosarcoma","Sarcoma")
  balanced_number <- 10000

} else if(layer_name == "layer_6_brain_nbm"){
  cells_keep <- c("Brain Cancer","Neuroblastoma")
  balanced_number <- 10000

}
#filter dataset
combined_metadata <- combined_metadata[which(combined_metadata$cell_type_merged_modified %in% cells_keep),]
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
for(j in index_for_loop){
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
  saveRDS(training, paste0("/unbalanced_set/unbalanced_training_count_mat_", layer_name,"_fold_", j, ".RDS"))
  saveRDS(training_metadata, paste0("/unbalanced_set/unbalanced_training_metadata_", layer_name,"_fold_", j, ".RDS"))
  saveRDS(testing, paste0("/unbalanced_set/unbalanced_testing_count_mat_", layer_name,"_fold_", j, ".RDS"))
  saveRDS(testing_metadata, paste0("/unbalanced_set/unbalanced_testing_metadata_", layer_name,"_fold_", j, ".RDS"))

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

  saveRDS(training, paste0("/split_by_folds/balanced_training_count_mat_", layer_name,"_fold_", j, ".RDS"))
  saveRDS(training_metadata, paste0("/split_by_folds/balanced_training_metadata_", layer_name,"_fold_", j, ".RDS"))
  print(j)
}






