#load dependencies
library(scATOMIC)
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)
library(cutoff)
library(copykat)
#load query dataset files
count_mat_file <- commandArgs(trailingOnly = "TRUE")
count_mat <- readRDS(paste0("testing_count_matrices/", count_mat_file))
metadata_name <- paste0(strsplit(count_mat_file, "count_mat.RDS")[[1]][1], "cell_metadata.RDS")
metadata_file <- readRDS(paste0("testing_metadata_files/", count_mat_file))
#filter low quality cells
pct_mt <- colSums(count_mat[grep("^MT-", row.names(count_mat)),])/colSums(count_mat) * 100
nFeatureRNA <- apply(count_mat, 2,function(x){length(which(x != 0))})
count_mat <- count_mat[, names(which(pct_mt < 25))]
count_mat <- count_mat[, intersect(names(which(nFeatureRNA > 500)), colnames(count_mat))]
#run scATOMIC
predictions <- run_scATOMIC(count_mat)
#run with CNV calculation for benchmarking figures and filtering
classifications <- create_summary_matrix(prediction_list = predictions,use_CNVs = T, modify_results = T, mc.cores = 20, raw_counts = count_matrix, min_prop = 0.5)
metadata <- metadata[row.names(classifications), ]
saveRDS(classifications, paste0("external_validation_results/", file_name, "_classification_with_CNV.RDS"))














