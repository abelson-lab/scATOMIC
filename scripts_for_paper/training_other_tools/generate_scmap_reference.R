#/usr/bin/Rscript
library(scmap)
library(dplyr)
library(data.table)
library(parallel)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
#set up chetah reference

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
cell_types <- gsub("CD4 CTL|CD4 Naive|CD4 Proliferating|CD4 TCM|CD4 TEM|Treg", "CD4+ T cell", cell_types)
cell_types <- gsub("CD8 Naive|CD8 Proliferating|CD8 TCM|CD8 TEM", "CD8+ T cell", cell_types)
cell_types <- gsub("NK|NK Proliferating|NK_CD56bright", "Natural killer cell", cell_types)
cell_types <- gsub("quiereferencent_stellate", "Pericytes and Stellate Cells", cell_types)

metadata_reference_matrix <- metadata_reference_matrix[-index_remove,]

metadata_reference_matrix$cell_type_merged <- cell_types


set.seed(123)
randomized_index <- sample(ncol(count_matrix_training))
randomized_dataset <- count_matrix_training[,randomized_index]
#10 fold cross validation
randomized_metadata <- metadata_reference_matrix[randomized_index,]
training <- randomized_dataset
training_metadata <- randomized_metadata
cell_types_merged <- levels(as.factor(as.character(training_metadata$cell_type_merged)))
colnames(training) <- row.names(training_metadata)
#split into 10 frames for or log normalization as the unsplit matrix is too big for R
training_list_index <- split(1:ncol(training), ceiling(seq_along(1:ncol(training))/(ncol(training)/20)))
training_list <- list()
for(i in 1:20){
  training_list[[i]] <- training[,training_list_index[[i]]]
  training_list[[i]] <- apply(training_list[[i]], 2, function(column) log2((column/sum(column) * 100000) + 1))
  training_list[[i]] <- as(as.matrix(training_list[[i]]), "sparseMatrix")

}
#merge split matrices
training_norm <- cbind(training_list[[1]], training_list[[2]], training_list[[3]], training_list[[4]],
                       training_list[[5]], training_list[[6]], training_list[[7]], training_list[[8]],
                       training_list[[9]], training_list[[10]],
                       training_list[[11]], training_list[[12]], training_list[[13]], training_list[[14]],
                       training_list[[15]], training_list[[16]], training_list[[17]], training_list[[18]],
                       training_list[[19]], training_list[[20]])


#create reference dataset format
reference <- SingleCellExperiment(assays = list(counts = training, logcounts = training_norm ),
                                  colData = DataFrame(celltypes = training_metadata$cell_type_merged))
#remove ribosomal genes
rowData(reference)$feature_symbol <- rownames(reference)
reference <- reference[!duplicated(rownames(reference)), ]
reference <- selectFeatures(reference, suppress_plot = TRUE)
rowData(reference)$scmap_features <- rowData(reference_old)$scmap_features
reference <- indexCell(reference)
saveRDS(reference, "reference_object_for_scmap.RDS")





