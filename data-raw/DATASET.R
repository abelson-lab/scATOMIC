library(data.table)
list_markers <- list.files("~/Downloads/expression_full/", full.names = T)
markers_master <- data.frame()
for(i in 1:length(list_markers)){
  markers_file <- data.table::fread(list_markers[i])
  markers_file <- as.data.frame(markers_file)
  cancer_type <- strsplit(list_markers[i], "diff_|.txt")[[1]][2]
  markers_file$cancer_type <- cancer_type
  markers_master <- rbind(markers_master, markers_file)
}





upreg_list <- markers_master[which(markers_master$`|log2FC|` > 1),]

downreg_list <- markers_master[which(markers_master$`|log2FC|` < -1),]


library(org.Hs.eg.db)
genes_rna_counts = upreg_list$Gene
annots = select(org.Hs.eg.db, keys = genes_rna_counts, columns=c( "ENSEMBL", "SYMBOL"), keytype = "ALIAS")
genes_rna_counts_ensemblID = c()
pb = txtProgressBar(min = 0, max = length(genes_rna_counts), style = 3)
for(i in 1:length(genes_rna_counts)){
  index_alias = which(annots$ALIAS == genes_rna_counts[i])
  index_gene_symbol = which(annots$SYMBOL == genes_rna_counts[i])

  if(length(index_alias)==1 & length(index_gene_symbol)==1){
    ensemblID = annots$ENSEMBL[index_alias]

  } else if(length(index_alias)>1 & length(index_gene_symbol)==1){
    ensemblID = annots$ENSEMBL[index_gene_symbol]

  } else if(length(index_gene_symbol) > 1){
    ensemblID = annots$ENSEMBL[index_gene_symbol[1]]
  } else if(length(index_gene_symbol) ==0 & length(index_alias) == 1){
    ensemblID = annots$ENSEMBL[index_alias]
  } else if(length(index_alias) > 1 & length(index_gene_symbol) == 0){
    ensemblID = annots$ENSEMBL[index_alias[1]]
  } else{
    print(paste0('check ',genes_rna_counts[i]))
  }
  if(!is.na(ensemblID)){
    genes_rna_counts_ensemblID[i] = ensemblID
  } else{
    genes_rna_counts_ensemblID[i] = NA
  }
  names(genes_rna_counts_ensemblID)[i] = genes_rna_counts[i]
  setTxtProgressBar(pb, i)
}
ensembl_IDs_rows = genes_rna_counts_ensemblID[genes_rna_counts]
index_na = which(is.na(ensembl_IDs_rows))
na_gene_names = names(genes_rna_counts_ensemblID)[index_na]
upreg_list = upreg_list[-which(upreg_list$Gene %in% na_gene_names),]

ensembl_IDs_rows = genes_rna_counts_ensemblID[upreg_list$Gene]
upreg_list$ensemblID = ensembl_IDs_rows



index_dup = which(duplicated(paste0(upreg_list$cancer_type, '___', upreg_list$ensemblID)))
#check for duplicates
if(length(index_dup)>0){
  upreg_list = upreg_list[-index_dup,]
}



genes_rna_counts = downreg_list$Gene
annots = select(org.Hs.eg.db, keys = genes_rna_counts, columns=c( "ENSEMBL", "SYMBOL"), keytype = "ALIAS")
genes_rna_counts_ensemblID = c()
pb = txtProgressBar(min = 0, max = length(genes_rna_counts), style = 3)
for(i in 1:length(genes_rna_counts)){
  index_alias = which(annots$ALIAS == genes_rna_counts[i])
  index_gene_symbol = which(annots$SYMBOL == genes_rna_counts[i])

  if(length(index_alias)==1 & length(index_gene_symbol)==1){
    ensemblID = annots$ENSEMBL[index_alias]

  } else if(length(index_alias)>1 & length(index_gene_symbol)==1){
    ensemblID = annots$ENSEMBL[index_gene_symbol]

  } else if(length(index_gene_symbol) > 1){
    ensemblID = annots$ENSEMBL[index_gene_symbol[1]]
  } else if(length(index_gene_symbol) ==0 & length(index_alias) == 1){
    ensemblID = annots$ENSEMBL[index_alias]
  } else if(length(index_alias) > 1 & length(index_gene_symbol) == 0){
    ensemblID = annots$ENSEMBL[index_alias[1]]
  } else{
    print(paste0('check ',genes_rna_counts[i]))
  }
  if(!is.na(ensemblID)){
    genes_rna_counts_ensemblID[i] = ensemblID
  } else{
    genes_rna_counts_ensemblID[i] = NA
  }
  names(genes_rna_counts_ensemblID)[i] = genes_rna_counts[i]
  setTxtProgressBar(pb, i)
}
ensembl_IDs_rows = genes_rna_counts_ensemblID[genes_rna_counts]
index_na = which(is.na(ensembl_IDs_rows))
na_gene_names = names(genes_rna_counts_ensemblID)[index_na]
downreg_list = downreg_list[-which(downreg_list$Gene %in% na_gene_names),]

ensembl_IDs_rows = genes_rna_counts_ensemblID[downreg_list$Gene]
downreg_list$ensemblID = ensembl_IDs_rows



index_dup = which(duplicated(paste0(downreg_list$cancer_type, '___', downreg_list$ensemblID)))
#check for ddownlicates
if(length(index_dup)>0){
  downreg_list = downreg_list[-index_dup,]
}











overlapping_genes <- intersect(upreg_list$ensemblID, downreg_list$ensemblID)
pan_cancer_upreg <- upreg_list[,"ensemblID"]
pan_cancer_upreg <- pan_cancer_upreg[-which(pan_cancer_upreg %in% overlapping_genes)]
pan_cancer_downreg <- downreg_list[,"ensemblID"]
pan_cancer_downreg <- pan_cancer_downreg[-which(pan_cancer_downreg %in% overlapping_genes)]








usethis::use_data(upreg_list, overwrite = TRUE)
usethis::use_data(downreg_list, overwrite = TRUE)
usethis::use_data(pan_cancer_upreg, overwrite = TRUE)
usethis::use_data(pan_cancer_downreg, overwrite = TRUE)


load(paste0("~/Documents/classifier_outputs_ensembl/","layer_1.RData"))
top_genes_unlisted_layer_1 <- top_genes_unlisted
model_layer_1 <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_1, overwrite = TRUE)
usethis::use_data(model_layer_1, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_2_blood.RData"))
top_genes_unlisted_layer_2_blood <- top_genes_unlisted
model_layer_2_blood <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_2_blood, overwrite = TRUE)
usethis::use_data(model_layer_2_blood, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_2_non_blood.RData"))
top_genes_unlisted_layer_2_normal_tissue_cancer <- top_genes_unlisted
model_layer_2_normal_tissue_cancer <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_2_normal_tissue_cancer, overwrite = TRUE)
usethis::use_data(model_layer_2_normal_tissue_cancer, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_3_B_cell.RData"))
top_genes_unlisted_layer_3_BCell <- top_genes_unlisted
model_layer_3_BCell <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_BCell, overwrite = TRUE)
usethis::use_data(model_layer_3_BCell, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_3_myeloid.RData"))
top_genes_unlisted_layer_3_MDC <- top_genes_unlisted
model_layer_3_MDC <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_MDC, overwrite = TRUE)
usethis::use_data(model_layer_3_MDC, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_3_non_stromal.RData"))

top_genes_unlisted_layer_3_non_stromal <- top_genes_unlisted
model_layer_3_non_stromal <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_non_stromal, overwrite = TRUE)
usethis::use_data(model_layer_3_non_stromal, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_3_stromal.RData"))
top_genes_unlisted_layer_3_stromal <- top_genes_unlisted
model_layer_3_stromal <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_stromal, overwrite = TRUE)
usethis::use_data(model_layer_3_stromal, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_3_TNK.RData"))
top_genes_unlisted_layer_3_TNK <- top_genes_unlisted
model_layer_3_TNK <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_TNK, overwrite = TRUE)
usethis::use_data(model_layer_3_TNK, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_4_CD4_CD8.RData"))
top_genes_unlisted_layer_4_CD4_CD8 <- top_genes_unlisted
model_layer_4_CD4_CD8 <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_CD4_CD8, overwrite = TRUE)
usethis::use_data(model_layer_4_CD4_CD8, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_4_CD8_NK.RData"))
top_genes_unlisted_layer_4_CD8_NK <- top_genes_unlisted
model_layer_4_CD8_NK <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_CD8_NK, overwrite = TRUE)
usethis::use_data(model_layer_4_CD8_NK, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_4_dendritic_cell.RData"))
top_genes_unlisted_layer_4_DC_cell <- top_genes_unlisted
model_layer_4_DC_cell <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_DC_cell, overwrite = TRUE)
usethis::use_data(model_layer_4_DC_cell, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_4_macrophage_monocyte.RData"))
top_genes_unlisted_layer_4_macrophage_monocyte <- top_genes_unlisted
model_layer_4_macrophage_monocyte <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_macrophage_monocyte, overwrite = TRUE)
usethis::use_data(model_layer_4_macrophage_monocyte, overwrite = TRUE)
rm(list = ls())


load(paste0("~/Documents/classifier_outputs_ensembl/","layer_4_GI.RData"))
top_genes_unlisted_layer_4_GI <- top_genes_unlisted
model_layer_4_GI <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_GI, overwrite = TRUE)
usethis::use_data(model_layer_4_GI, overwrite = TRUE)
rm(list = ls())
####
load(paste0("~/Documents/classifier_outputs_ensembl/","layer_4_non_GI.RData"))
top_genes_unlisted_layer_4_non_GI <- top_genes_unlisted
model_layer_4_non_GI <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_non_GI, overwrite = TRUE)
usethis::use_data(model_layer_4_non_GI, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_4_soft_tissue_neuro.RData"))
top_genes_unlisted_layer_4_soft_neuro_cancer <- top_genes_unlisted
model_layer_4_soft_neuro_cancer <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_soft_neuro_cancer, overwrite = TRUE)
usethis::use_data(model_layer_4_soft_neuro_cancer, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_biliary.RData"))
top_genes_unlisted_layer_5_biliary <- top_genes_unlisted
model_layer_5_biliary <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_biliary, overwrite = TRUE)
usethis::use_data(model_layer_5_biliary, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_breast_lung_prostate.RData"))
top_genes_unlisted_layer_5_breast_lung_prostate <- top_genes_unlisted
model_layer_5_breast_lung_prostate <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_breast_lung_prostate, overwrite = TRUE)
usethis::use_data(model_layer_5_breast_lung_prostate, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_digestive.RData"))
top_genes_unlisted_layer_5_digestive <- top_genes_unlisted
model_layer_5_digestive <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_digestive, overwrite = TRUE)
usethis::use_data(model_layer_5_digestive, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_ov_endo_kid.RData"))
top_genes_unlisted_layer_5_ov_endo_kid <- top_genes_unlisted
model_layer_5_ov_endo_kid  <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_ov_endo_kid, overwrite = TRUE)
usethis::use_data(model_layer_5_ov_endo_kid, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_soft_tissue_neuro.RData"))
top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin <- top_genes_unlisted
model_layer_5_soft_neuro_cancer_no_lung_skin <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin, overwrite = TRUE)
usethis::use_data(model_layer_5_soft_neuro_cancer_no_lung_skin, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_6_brain_nbm.RData"))
top_genes_unlisted_layer_6_brain_nbm <- top_genes_unlisted
model_layer_6_brain_nbm <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_6_brain_nbm, overwrite = TRUE)
usethis::use_data(model_layer_6_brain_nbm, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_6_soft_tissue.RData"))
top_genes_unlisted_layer_6_soft_tissue <- top_genes_unlisted

model_layer_6_soft_tissue <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_6_soft_tissue, overwrite = TRUE)
usethis::use_data(model_layer_6_soft_tissue, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_6_breast.RData"))
top_genes_unlisted_layer_6_breast <- top_genes_unlisted

model_layer_6_breast <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_6_breast, overwrite = TRUE)
usethis::use_data(model_layer_6_breast, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_CD4.RData"))
top_genes_unlisted_layer_5_CD4 <- top_genes_unlisted
model_layer_5_CD4 <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_CD4, overwrite = TRUE)
usethis::use_data(model_layer_5_CD4, overwrite = TRUE)
rm(list = ls())


load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_CD8.RData"))
top_genes_unlisted_layer_5_CD8 <- top_genes_unlisted
model_layer_5_CD8 <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_CD8, overwrite = TRUE)
usethis::use_data(model_layer_5_CD8, overwrite = TRUE)
rm(list = ls())

demo_lung_data <- demo_lung_data

usethis::use_data(demo_lung_data, overwrite = TRUE)
rm(list = ls())






load(paste0("~/Documents/classifier_outputs_ensembl/","layer_4_CAFs.RData"))
top_genes_unlisted_layer_4_CAFs <- top_genes_unlisted
model_layer_4_CAFs <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_4_CAFs, overwrite = TRUE)
usethis::use_data(model_layer_4_CAFs, overwrite = TRUE)
rm(list = ls())


load(paste0("~/Documents/classifier_outputs_ensembl","/layer_5_cDC.RData"))
top_genes_unlisted_layer_5_cDC <- top_genes_unlisted
model_layer_5_cDC <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_cDC, overwrite = TRUE)
usethis::use_data(model_layer_5_cDC, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_macrophage.RData"))
top_genes_unlisted_layer_6_macrophage <- top_genes_unlisted
model_layer_6_macrophage <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_6_macrophage, overwrite = TRUE)
usethis::use_data(model_layer_6_macrophage, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/classifier_outputs_ensembl/","layer_5_monocyte.RData"))
top_genes_unlisted_layer_6_monocyte <- top_genes_unlisted
model_layer_6_monocyte <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_6_monocyte, overwrite = TRUE)
usethis::use_data(model_layer_6_monocyte, overwrite = TRUE)
rm(list = ls())


load(paste0("~/Documents/neutrophil_classifier/outs/ensembl_classifier_layer_5_macrophage_neutrophil.RData"))
top_genes_unlisted_layer_5_macrophage_neutrophil <- top_genes_unlisted
model_layer_5_macrophage_neutrophil <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_macrophage_neutrophil, overwrite = TRUE)
usethis::use_data(model_layer_5_macrophage_neutrophil, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/neutrophil_classifier/outs/ensembl_classifier_layer_5_monocyte_neutrophil.RData"))
top_genes_unlisted_layer_5_monocyte_neutrophil <- top_genes_unlisted
model_layer_5_monocyte_neutrophil <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_monocyte_neutrophil, overwrite = TRUE)
usethis::use_data(model_layer_5_monocyte_neutrophil, overwrite = TRUE)
rm(list = ls())




