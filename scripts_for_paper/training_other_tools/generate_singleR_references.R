#get reference for singleR
#/usr/bin/Rscript
library(SingleR)
library(parallel)
library(Matrix)
library(CHETAH)
library(dplyr)
library(data.table)
library(Seurat)
library(SingleCellExperiment)
library(scater)
count_matrix_training <- readRDS("combined_RNA_count_matrix_layer_1.RDS")
metadata_reference_matrix <- readRDS("combined_RNA_count_matrix_layer_1_metadata.RDS")
cell_types <- metadata_reference_matrix[,"cell_type_merged"]
index_remove <- which(cell_types %in%c("Doublet", "dnT", "ILC","CD4 CTL", "gdT","MAIT","Proliferating" , "Head and Neck Cancer", "Thyroid Cancer",
                                       "Eryth", "Platelet"))
count_matrix_training <- count_matrix_training[,-index_remove]
cell_types <- cell_types[-index_remove]
cell_types <- gsub("CD14 Mono|CD16 Mono", "Macrophage", cell_types)
cell_types <- gsub("ASDC|cDC1|cDC2", "Dendritic Cell", cell_types)

cell_types <- gsub("B intermediate|B memory|B naive", "B cell", cell_types)
cell_types <- gsub("CD4 Naive|CD4 Proliferating|CD4 TCM|CD4 TEM|Treg", "CD4+ T cell", cell_types)
cell_types <- gsub("CD8 Naive|CD8 Proliferating|CD8 TCM|CD8 TEM", "CD8+ T cell", cell_types)
cell_types <- gsub("NK|NK Proliferating|NK_CD56bright", "Natural killer cell", cell_types)




metadata_reference_matrix <- metadata_reference_matrix[-index_remove,]
#need to load the new gene set

metadata_reference_matrix$cell_type_merged <- cell_types
#sample 2500 of each
set.seed(123)
randomized_index <- sample(ncol(count_matrix_training))
randomized_dataset <- count_matrix_training[,randomized_index]
#10 fold cross validation
randomized_metadata <- metadata_reference_matrix[randomized_index,]



training <- randomized_dataset
training_metadata <- randomized_metadata
cell_types_merged <- levels(as.factor(as.character(training_metadata$cell_type_merged)))

colnames(training) <- row.names(training_metadata)


print("getting sceref")
sceRef <- SingleCellExperiment(list(counts=training), colData=DataFrame(label=training_metadata$cell_type_merged))


print("getting logcounts")

sceRef <- logNormCounts(sceRef)
print("saving")
saveRDS(sceRef, "singleR_reference.RDS")
#part 2 generate dataset specific reference reflective of overlapping features in query dataset and reference
#/usr/bin/Rscript
library(SingleR)
library(parallel)
library(Matrix)
library(dplyr)
library(data.table)
library(Seurat)
library(SingleCellExperiment)
library(scater)
sceRef <- readRDS("singleR_reference.RDS")
#build singleR pretrained for each dataset
dataset_index <- as.numeric(commandArgs(trailingOnly = T)[1])
#dataset_index <- 1
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

list_all_samples <- list.files("split_by_patient_matrices/", pattern = dataset_use, full.names = T)
for(i in 1:length(list_all_samples)){
  mat <- readRDS(list_all_samples[i])
  if(i ==1 ){
    master_gene_list <-row.names(mat)
  } else{
    master_gene_list <- intersect(master_gene_list, row.names(mat))
  }
}

common <- intersect(rownames(sceRef), master_gene_list)
print("training aggregated")
set.seed(2000)
trained_agg <- trainSingleR(sceRef[common,], labels=sceRef$label, aggr.ref=TRUE)
saveRDS(trained_agg, paste0("singleR_references_aggregated_ref/", dataset_use, "_singleR_trained_reference_object.RDS"))


