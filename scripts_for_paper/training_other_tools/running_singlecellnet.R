#/usr/bin/Rscript
library(Matrix)
library(dplyr)
library(data.table)
library(Seurat)
library(singleCellNet)
#pass on biopsy file as argument
file_name <- commandArgs(trailingOnly = T)[1]
#file name starts with dataset name
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
print(file_name)
dataset_use <- c()
for(i in 1:length(datasets)){
  if(grepl(datasets[i], file_name)){
    print(datasets[i])
    dataset_use <- datasets[i]
  }
}
#load dataset specific Singlecellnet model
reference_file <- list.files("single_cell_net_references", full.names = T, pattern = dataset_use)
SCN_reference <- readRDS(reference_file)

#load query file
metadata_name <- paste0(strsplit(file_name, "count_mat.RDS")[[1]][1], "cell_metadata.RDS")
path_to_matrix <- paste0("split_by_patient_matrices/", file_name)
path_to_metadata <- paste0("split_by_patient_metadata/", metadata_name)
expQuery <- readRDS(path_to_matrix)
metadata <- readRDS(path_to_metadata)
#filter low quality cells
pct_mt <- colSums(expQuery[grep("^MT-", row.names(expQuery)),])/colSums(expQuery) * 100
nFeatureRNA <- apply(expQuery, 2,function(x){length(which(x != 0))})
expQuery <- expQuery[, names(which(pct_mt < 25))]
expQuery <- expQuery[, intersect(names(which(nFeatureRNA > 500)), colnames(expQuery))]
metadata <- metadata[colnames(expQuery),]
#load ref mat to find intersecting genes
count_matrix_training <- readRDS("/.mounts/labs/abelsonlab/public/Ido/july_2021_retrain_for_paper/training_datasets/combined_RNA_count_matrix_layer_1_july15.RDS")

commonGenes<-intersect(rownames(count_matrix_training), rownames(expQuery))
expQuery <- expQuery[commonGenes, ]
rm(count_matrix_training)



stQuery <- metadata[colnames(expQuery),]
#run prediction
system.time(classsification <- scn_predict(SCN_reference[['cnProc']], expQuery, nrand = 50))

classsification <- classsification[,row.names(stQuery)]
stQuery <- assign_cate(classRes = classsification, sampTab = stQuery, cThresh = 0.5)
metadata$SCN_classification <- stQuery[row.names(metadata), "category"]

saveRDS(metadata, paste0("single_cell_net_results/", file_name, "_results_singlecellnet.RDS"))















