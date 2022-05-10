#/usr/bin/Rscript
#run scATOMIC on test sets
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
library(scATOMIC)
library(copykat)
#run from the directory
#load testing set file from script 1
file_name <- commandArgs(trailingOnly = T)[1]
fold  <- commandArgs(trailingOnly = T)[2]
metadata_name <- paste0("unbalanced_testing_metadata_", strsplit(file_name, "unbalanced_testing_count_mat_")[[1]][2])
print(file_name)
path_to_matrix <- paste0("testing_sets/", file_name)
path_to_metadata <- paste0("testing_sets/", metadata_name)
count_matrix <- readRDS(path_to_matrix)
metadata <- readRDS(path_to_metadata)
#QC filter testing set
pct_mt <- colSums(count_matrix[grep("^MT-", row.names(count_matrix)),])/colSums(count_matrix) * 100
nFeatureRNA <- apply(count_matrix, 2,function(x){length(which(x != 0))})
count_matrix <- count_matrix[, names(which(pct_mt < 25))]
count_matrix <- count_matrix[, intersect(names(which(nFeatureRNA > 500)), colnames(count_matrix))]

#load k fold models
load(paste0("classifier_outputs_unbalanced_per_fold/layer_1_fold_",fold,".RData"))
top_genes_unlisted_layer_1 <- top_genes_unlisted
model_layer_1 <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_2_blood_fold_",fold,".RData"))
top_genes_unlisted_layer_2_blood <- top_genes_unlisted
model_layer_2_blood <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_2_non_blood_fold_",fold,".RData"))
top_genes_unlisted_layer_2_normal_tissue_cancer <- top_genes_unlisted
model_layer_2_normal_tissue_cancer <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_3_B_cell_fold_",fold,".RData"))
top_genes_unlisted_layer_3_BCell <- top_genes_unlisted
model_layer_3_BCell <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_3_MDC_fold_",fold,".RData"))
top_genes_unlisted_layer_3_MDC <- top_genes_unlisted
model_layer_3_MDC <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_3_non_stromal_fold_",fold,".RData"))
top_genes_unlisted_layer_3_non_stromal <- top_genes_unlisted
model_layer_3_non_stromal <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_3_stromal_fold_",fold,".RData"))
top_genes_unlisted_layer_3_stromal <- top_genes_unlisted
model_layer_3_stromal <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_3_TNK_fold_",fold,".RData"))
top_genes_unlisted_layer_3_TNK <- top_genes_unlisted
model_layer_3_TNK <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_4_CD4_CD8_fold_",fold,".RData"))
top_genes_unlisted_layer_4_CD4_CD8 <- top_genes_unlisted
model_layer_4_CD4_CD8 <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_4_CD8_NK_cell_fold_",fold,".RData"))
top_genes_unlisted_layer_4_CD8_NK <- top_genes_unlisted
model_layer_4_CD8_NK <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_4_DC_fold_",fold,".RData"))
top_genes_unlisted_layer_4_DC_cell <- top_genes_unlisted
model_layer_4_DC_cell <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_4_GI_fold_",fold,".RData"))
top_genes_unlisted_layer_4_GI <- top_genes_unlisted
model_layer_4_GI <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_4_non_GI_fold_",fold,".RData"))
top_genes_unlisted_layer_4_non_GI <- top_genes_unlisted
model_layer_4_non_GI <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_4_soft_tissue_neuro_fold_",fold,".RData"))
top_genes_unlisted_layer_4_soft_neuro_cancer <- top_genes_unlisted
model_layer_4_soft_neuro_cancer <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_5_biliary_fold_",fold,".RData"))
top_genes_unlisted_layer_5_biliary <- top_genes_unlisted
model_layer_5_biliary <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_5_breast_lung_prostate_fold_",fold,".RData"))
top_genes_unlisted_layer_5_breast_lung_prostate <- top_genes_unlisted
model_layer_5_breast_lung_prostate <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_5_digestive_fold_",fold,".RData"))
top_genes_unlisted_layer_5_digestive <- top_genes_unlisted
model_layer_5_digestive <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_5_ov_endo_kid_fold_",fold,".RData"))
top_genes_unlisted_layer_5_ov_endo_kid <- top_genes_unlisted
model_layer_5_ov_endo_kid  <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_5_soft_tissue_neuro_fold_",fold,".RData"))
top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin <- top_genes_unlisted
model_layer_5_soft_neuro_cancer_no_lung_skin <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_6_brain_nbm_fold_",fold,".RData"))
top_genes_unlisted_layer_6_brain_nbm <- top_genes_unlisted
model_layer_6_brain_nbm <- rf_classifier_cell_lines

load(paste0("classifier_outputs_unbalanced_per_fold/layer_6_soft_tissue_fold_",fold,".RData"))
top_genes_unlisted_layer_6_soft_tissue <- top_genes_unlisted
model_layer_6_soft_tissue <- rf_classifier_cell_lines



#run scATOMIC without cutoffs
predictions <- run_scATOMIC(rna_counts = count_matrix,confidence_cutoff = F)

saveRDS(predictions, paste0("results_testing_sets/", file_name, "_predictions_forced.RDS"))
#create summary matrix
classifications <- create_summary_matrix(prediction_list = predictions,modify_results = F,raw_counts = count_matrix, confidence_cutoff = F)
metadata <- metadata[row.names(classifications), ]
classifications <- cbind(classifications, metadata)
saveRDS(classifications, paste0("results_testing_sets/", file_name, "_classification_no_cutoff.RDS"))




