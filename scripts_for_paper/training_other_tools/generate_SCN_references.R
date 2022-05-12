#/usr/bin/Rscript
#training singlecellnet
library(Matrix)
library(dplyr)
library(data.table)
library(Seurat)
library(singleCellNet)
#load reference
count_matrix_training <- readRDS("combined_RNA_count_matrix_layer_1.RDS")
metadata_reference_matrix <- readRDS("/.mounts/labs/abelsonlab/public/Ido/july_2021_retrain_for_paper/training_datasets/combined_RNA_count_matrix_layer_1_july15_metadata.RDS")
cell_types <- metadata_reference_matrix[,"cell_type_merged"]
index_remove <- which(cell_types %in%c("Doublet", "dnT", "ILC","CD4 CTL", "gdT","MAIT","Proliferating" , "Head and Neck Cancer", "Thyroid Cancer",
                                       "Eryth", "Platelet"))
count_matrix_training <- count_matrix_training[,-index_remove]
cell_types <- cell_types[-index_remove]
cell_types <- gsub("CD14 Mono|CD16 Mono", "Macrophage", cell_types)
cell_types <- gsub("ASDC|cDC1|cDC2", "Dendritic Cell", cell_types)
cell_types <- gsub("B intermediate|B memory|B naive", "B cell", cell_types)
cell_types <- gsub("CD4 CTL|CD4 Naive|CD4 Proliferating|CD4 TCM|CD4 TEM|Treg", "CD4+ T cell", cell_types)
cell_types <- gsub("CD8 Naive|CD8 Proliferating|CD8 TCM|CD8 TEM", "CD8+ T cell", cell_types)
cell_types <- gsub("NK|NK Proliferating|NK_CD56bright", "Natural killer cell", cell_types)
cell_types <- gsub("quiescent_stellate", "Pericytes and Stellate Cells", cell_types)

metadata_reference_matrix <- metadata_reference_matrix[-index_remove,]

metadata_reference_matrix$cell_type_merged <- cell_types




#create training matrix in SCN format
stTM <- metadata_reference_matrix
expTMraw <- count_matrix_training

stList<-splitCommon(sampTab = stTM, ncells = 2500, dLevel = "cell_type_merged" )

stTrain<-stList[[1]]
stTrain$cell <- row.names(stTrain)
expTrain <- expTMraw[,as.character(row.names(stTrain))]




#trained a SCN model for each datatset
#each query dataset has a slightly different gene list so we trained a model reflecting the overlapping refrence and query genes for each dataset
dataset_index <- as.numeric(commandArgs(trailingOnly = T)[1])
datasets <- c("Bi_et_al_kidney",
              "Chen_et_al_bladder",
              "Chen_et_al_prostate",
              "Couterier_et_al_GBM",
              "Dong_et_al_nbm",
              "Gao_et_al_copykat",
              "Jerby-Arnon_et_al_melanoma",
              "Kim_et_al_lung",
              "Lee_et_al_CRC",
              "Ma_et_al_liver",
              "Peng_et_al_pancreas",
              "Qian_et_al_pan_cancer",
              "Slyper_et_al_pan_cancer",
              "Wu_et_al_breast_atlas",
              "Wu_et_al_pan_cancer",
              "Young_et_al_kidney")
dataset_use <- datasets[dataset_index]
print(dataset_use)
file_name <- list.files("split_by_patient_matrices/", pattern = dataset_use)[1]


print(file_name)
path_to_matrix <- paste0("split_by_patient_matrices/", file_name)
expQuery <- readRDS(path_to_matrix)
list_all_samples <- list.files("split_by_patient_matrices/", pattern = dataset_use, full.names = T)
for(i in 1:length(list_all_samples)){
  mat <- readRDS(list_all_samples[i])
  if(i ==1 ){
    master_gene_list <-row.names(mat)
  } else{
    master_gene_list <- intersect(master_gene_list, row.names(mat))
  }
}
#retrieve gene list for SCN model
commonGenes<-intersect(rownames(expTrain), master_gene_list)
expTrain <- expTrain[commonGenes, ]
expQuery <- expQuery[commonGenes, ]
#train SCN model
print("Training")
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "cell_type_merged", colName_samp = "cell"))
print("Training Done")
#save SCN reference model
saveRDS(class_info, paste0("/single_cell_net_references/", dataset_use, "_reference_object.RDS"))
