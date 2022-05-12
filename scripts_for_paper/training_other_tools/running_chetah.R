#/usr/bin/Rscript
library(parallel)
library(Matrix)
library(CHETAH)
library(dplyr)
library(data.table)
library(Seurat)
library(SingleCellExperiment)
#pass on biopsy file as argument
file_name <- commandArgs(trailingOnly = T)[1]
print(file_name)
#load query file and metadata
metadata_name <- paste0(strsplit(file_name, "count_mat.RDS")[[1]][1], "cell_metadata.RDS")
path_to_matrix <- paste0("split_by_patient_matrices/", file_name)
path_to_metadata <- paste0("split_by_patient_metadata/", metadata_name)
count_matrix <- readRDS(path_to_matrix)
metadata <- readRDS(path_to_metadata)
#filter low quality cells
pct_mt <- colSums(count_matrix[grep("^MT-", row.names(count_matrix)),])/colSums(count_matrix) * 100
nFeatureRNA <- apply(count_matrix, 2,function(x){length(which(x != 0))})
count_matrix <- count_matrix[, names(which(pct_mt < 25))]
count_matrix <- count_matrix[, intersect(names(which(nFeatureRNA > 500)), colnames(count_matrix))]
#load chetah reference
reference_chetah <- readRDS("reference_object_for_chetah.RDS")
input_chetah <- SingleCellExperiment(assays = list(counts = count_matrix))
print("Starting CHETAH")
input_chetah <- CHETAHclassifier(input = input_chetah, ref_cells = reference_chetah, ref_c = "logcounts",thresh = 0)
chetah_annotations <- input_chetah@colData$celltype_CHETAH
print("Done CHETAH")
metadata <- metadata[names(chetah_annotations), ]

metadata <- cbind(metadata, chetah_annotations)

saveRDS(metadata, paste0("chetah_results/", file_name, "_results_chetah.RDS"))









