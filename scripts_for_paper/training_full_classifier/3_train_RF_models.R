#/usr/bin/Rscript
#train the random forest models for each branch
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(dplyr)
library(Seurat)
library(Matrix)
#pass the balanced count matrix generated in script 1 as the arg
file_name <- commandArgs(trailingOnly = T)[1]
#set cutoff off for filtering DGE list, default is 0 and 0.4
number_sd <- 0
percent_cutoff <- 0.4
metadata_name <- paste0("balanced_training_metadata_",strsplit(file_name, "balanced_training_count_mat_")[[1]][2])
#load files
training <- readRDS(paste0("balanced_set_full_classifier/", file_name))
training_metadata <- readRDS(paste0("balanced_set_full_classifier/", metadata_name))
marker_file_pattern <- strsplit(file_name, "balanced_training_count_mat_|.RDS")[[1]][2]
file_list <- list.files("markers_full_classifier/", pattern = marker_file_pattern, full.names = T)
per_dataset_marker_list <- readRDS(file_list)
#filter DGEs to selected features
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
    top_pos_per_test[[m]] <- per_dataset_marker_list[[m]] %>% top_n(n = 50, wt = diff)
  }
  if(nrow(top_pos_per_test[[m]]) > 200){
    top_pos_per_test[[m]] <- top_pos_per_test[[m]] %>% top_n(n = 200, wt = diff)
  }


  top_pos_per_test[[m]] <- top_pos_per_test[[m]]$gene
}
names(top_pos_per_test)  <- names(per_dataset_marker_list)
top_pos_per_test <- as.vector(unlist(top_pos_per_test, use.names = FALSE))
top_genes_unlisted <- unique(top_pos_per_test)
# library size normalize training matrix
training <- t(training)
training <- library.size.normalize(training)
training <-t(training)
#filter matrix to include only selected features
training <- na.omit(training[top_genes_unlisted,])
#transform matrix from reads to relative reads
training <- as.data.frame(apply(training, 2, function(x){
  x/sum(x)
}))
#format training
training <- as.data.frame(t(training))
cell_names <- row.names(training)
training_metadata <- training_metadata[cell_names,]
cell_class <- as.factor(training_metadata$cell_type_merged_modified)
training <- cbind(cell_class, training)
names_genes <- colnames(training)
names_genes <- gsub("-", "_", names_genes)
colnames(training) <- names_genes
training <- na.omit(training)
#train RF model, using 5 cores here - change to cores to desired
print("training started")
library(foreach)
library(doParallel)
registerDoParallel(cores=5)
#trained forest with 5 * 100 trees
rf_classifier_cell_lines <- foreach(ntree=rep(100, 5), .combine=randomForest::combine,
                                    .multicombine=TRUE, .packages='randomForest') %dopar% {
                                      randomForest(cell_class ~ ., data=training, ntree=ntree)
                                    }
print("training ended")
# keep genes expressed in at least 10 cells

#save model
save(list = c("rf_classifier_cell_lines", "top_genes_unlisted" ), file = paste0("classifier_outputs_full_classifier/", marker_file_pattern, ".RData"))

#trained models and gene lists are included in this R package in the data folder



