#/usr/bin/Rscript
library(parallel)
library(Matrix)
library(scmap)
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
#load scmap reference
scmap_reference <- readRDS("reference_object_for_scmap.RDS")
#normalize query matrix
count_matrix_norm  <- apply(count_matrix, 2, function(column) log2((column/sum(column) * 100000) + 1))
count_matrix_norm <- as(as.matrix(count_matrix_norm), "sparseMatrix")
#create sce object for query
query_sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=count_matrix, logcounts = count_matrix_norm))
# add feature_symbol column (i.e. the gene symbols)
rowData(query_sce)$feature_symbol <- rownames(query_sce)
scmap_cell_reference <- metadata(scmap_reference)$scmap_cell_index
#run scmap-cell
print("Starting scmap-cell")

nearest_neighbours <- scmap::scmapCell(projection=query_sce,
                                       index_list = list(ref = scmap_cell_reference))
print("done scmap-cell")
#add metadata
scmap_cell_metadata <- colData(scmap_reference)
colnames(scmap_cell_metadata) <- "celltypes"
mode_label <- function(neighbours, metadata=scmap_cell_metadata$celltypes) {
  freq <- table(metadata[neighbours])
  label <- names(freq)[which(freq == max(freq))]
  if (length(label) > 1) {return("Unassigned")}
  return(label)
}
# Apply these labels to the query cells
scmap_cell_labs <- apply(nearest_neighbours$ref$cells, 2, mode_label)
# Add the labels to the query object
colData(query_sce)$scmap_cell <- scmap_cell_labs

scmap_annotations <- query_sce@colData$scmap_cell

metadata <- metadata[names(scmap_annotations), ]

metadata <- cbind(metadata, scmap_annotations)
saveRDS(metadata, paste0("scmap_cell_results/", file_name, "_results_scmap_cell.RDS"))




