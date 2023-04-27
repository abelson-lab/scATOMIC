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

overlapping_genes <- intersect(upreg_list$Gene, downreg_list$Gene)
pan_cancer_upreg <- upreg_list[,"Gene"]
pan_cancer_upreg <- pan_cancer_upreg[-which(pan_cancer_upreg %in% overlapping_genes)]
pan_cancer_downreg <- downreg_list[,"Gene"]
pan_cancer_downreg <- pan_cancer_downreg[-which(pan_cancer_downreg %in% overlapping_genes)]
usethis::use_data(upreg_list, overwrite = TRUE)
usethis::use_data(downreg_list, overwrite = TRUE)
usethis::use_data(pan_cancer_upreg, overwrite = TRUE)
usethis::use_data(pan_cancer_downreg, overwrite = TRUE)


load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_1.RData"))
top_genes_unlisted_layer_1 <- top_genes_unlisted
model_layer_1 <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_1, overwrite = TRUE)
usethis::use_data(model_layer_1, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_2_blood.RData"))
top_genes_unlisted_layer_2_blood <- top_genes_unlisted
model_layer_2_blood <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_2_blood, overwrite = TRUE)
usethis::use_data(model_layer_2_blood, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_2_non_blood.RData"))
top_genes_unlisted_layer_2_normal_tissue_cancer <- top_genes_unlisted
model_layer_2_normal_tissue_cancer <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_2_normal_tissue_cancer, overwrite = TRUE)
usethis::use_data(model_layer_2_normal_tissue_cancer, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_3_B_cell.RData"))
top_genes_unlisted_layer_3_BCell <- top_genes_unlisted
model_layer_3_BCell <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_BCell, overwrite = TRUE)
usethis::use_data(model_layer_3_BCell, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_3_myeloid.RData"))
top_genes_unlisted_layer_3_MDC <- top_genes_unlisted
model_layer_3_MDC <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_MDC, overwrite = TRUE)
usethis::use_data(model_layer_3_MDC, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_3_non_stromal <- scATOMIC::top_genes_unlisted_layer_3_non_stromal
model_layer_3_non_stromal <- scATOMIC::model_layer_3_non_stromal

usethis::use_data(top_genes_unlisted_layer_3_non_stromal, overwrite = TRUE)
usethis::use_data(model_layer_3_non_stromal, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_3_stromal.RData"))
top_genes_unlisted_layer_3_stromal <- top_genes_unlisted
model_layer_3_stromal <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_stromal, overwrite = TRUE)
usethis::use_data(model_layer_3_stromal, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_3_TNK.RData"))
top_genes_unlisted_layer_3_TNK <- top_genes_unlisted
model_layer_3_TNK <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_TNK, overwrite = TRUE)
usethis::use_data(model_layer_3_TNK, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_4_CD4_CD8.RData"))
top_genes_unlisted_layer_4_CD4_CD8 <- top_genes_unlisted
model_layer_4_CD4_CD8 <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_CD4_CD8, overwrite = TRUE)
usethis::use_data(model_layer_4_CD4_CD8, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_4_CD8_NK.RData"))
top_genes_unlisted_layer_4_CD8_NK <- top_genes_unlisted
model_layer_4_CD8_NK <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_CD8_NK, overwrite = TRUE)
usethis::use_data(model_layer_4_CD8_NK, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_4_dendritic_cell.RData"))
top_genes_unlisted_layer_4_DC_cell <- top_genes_unlisted
model_layer_4_DC_cell <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_DC_cell, overwrite = TRUE)
usethis::use_data(model_layer_4_DC_cell, overwrite = TRUE)
rm(list = ls())
load(paste0("~/Documents/scATOMIC_revisions_with_PBMC/classifier_outputs_scATOMIC_revisions_with_PBMC_v2/", "Chen_et_al_33420488", "_heldout_layer_4_macrophage_monocyte.RData"))
top_genes_unlisted_layer_4_macrophage_monocyte <- top_genes_unlisted
model_layer_4_macrophage_monocyte <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_macrophage_monocyte, overwrite = TRUE)
usethis::use_data(model_layer_4_macrophage_monocyte, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_4_GI <- scATOMIC::top_genes_unlisted_layer_4_GI
model_layer_4_GI <- scATOMIC::model_layer_4_GI

usethis::use_data(top_genes_unlisted_layer_4_GI, overwrite = TRUE)
usethis::use_data(model_layer_4_GI, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_4_non_GI <- scATOMIC::top_genes_unlisted_layer_4_non_GI
model_layer_4_non_GI <- scATOMIC::model_layer_4_non_GI

usethis::use_data(top_genes_unlisted_layer_4_non_GI, overwrite = TRUE)
usethis::use_data(model_layer_4_non_GI, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_4_soft_neuro_cancer <- scATOMIC::top_genes_unlisted_layer_4_soft_neuro_cancer
model_layer_4_soft_neuro_cancer <- scATOMIC::model_layer_4_soft_neuro_cancer

usethis::use_data(top_genes_unlisted_layer_4_soft_neuro_cancer, overwrite = TRUE)
usethis::use_data(model_layer_4_soft_neuro_cancer, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_5_biliary <- scATOMIC::top_genes_unlisted_layer_5_biliary
model_layer_5_biliary <- scATOMIC::model_layer_5_biliary

usethis::use_data(top_genes_unlisted_layer_5_biliary, overwrite = TRUE)
usethis::use_data(model_layer_5_biliary, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_5_breast_lung_prostate <- scATOMIC::top_genes_unlisted_layer_5_breast_lung_prostate
model_layer_5_breast_lung_prostate <- scATOMIC::model_layer_5_breast_lung_prostate

usethis::use_data(top_genes_unlisted_layer_5_breast_lung_prostate, overwrite = TRUE)
usethis::use_data(model_layer_5_breast_lung_prostate, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_5_digestive <- scATOMIC::top_genes_unlisted_layer_5_digestive
model_layer_5_digestive <- scATOMIC::model_layer_5_digestive

usethis::use_data(top_genes_unlisted_layer_5_digestive, overwrite = TRUE)
usethis::use_data(model_layer_5_digestive, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_5_ov_endo_kid <- scATOMIC::top_genes_unlisted_layer_5_ov_endo_kid
model_layer_5_ov_endo_kid  <- scATOMIC::model_layer_5_ov_endo_kid

usethis::use_data(top_genes_unlisted_layer_5_ov_endo_kid, overwrite = TRUE)
usethis::use_data(model_layer_5_ov_endo_kid, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin <- scATOMIC::top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin
model_layer_5_soft_neuro_cancer_no_lung_skin <- scATOMIC::model_layer_5_soft_neuro_cancer_no_lung_skin

usethis::use_data(top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin, overwrite = TRUE)
usethis::use_data(model_layer_5_soft_neuro_cancer_no_lung_skin, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_6_brain_nbm <- scATOMIC::top_genes_unlisted_layer_6_brain_nbm
model_layer_6_brain_nbm <- scATOMIC::model_layer_6_brain_nbm
usethis::use_data(top_genes_unlisted_layer_6_brain_nbm, overwrite = TRUE)
usethis::use_data(model_layer_6_brain_nbm, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_6_soft_tissue <- scATOMIC::top_genes_unlisted_layer_6_soft_tissue

model_layer_6_soft_tissue <- scATOMIC::model_layer_6_soft_tissue
usethis::use_data(top_genes_unlisted_layer_6_soft_tissue, overwrite = TRUE)
usethis::use_data(model_layer_6_soft_tissue, overwrite = TRUE)
rm(list = ls())
top_genes_unlisted_layer_6_breast <- scATOMIC::top_genes_unlisted_layer_6_breast

model_layer_6_breast <- scATOMIC::model_layer_6_breast
usethis::use_data(top_genes_unlisted_layer_6_breast, overwrite = TRUE)
usethis::use_data(model_layer_6_breast, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Downloads/data_expression_T_cell_atlas/classifier_outputs_T_cells/classifier_layer_5_CD4_subtypes.RData"))
top_genes_unlisted_layer_5_CD4 <- top_genes_unlisted
model_layer_5_CD4 <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_CD4, overwrite = TRUE)
usethis::use_data(model_layer_5_CD4, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Downloads/data_expression_T_cell_atlas/classifier_outputs_T_cells/classifier_layer_5_CD8_subtypes.RData"))
top_genes_unlisted_layer_5_CD8 <- top_genes_unlisted
model_layer_5_CD8 <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_CD8, overwrite = TRUE)
usethis::use_data(model_layer_5_CD8, overwrite = TRUE)
rm(list = ls())

demo_lung_data <- scATOMIC::demo_lung_data

usethis::use_data(demo_lung_data, overwrite = TRUE)
rm(list = ls())






load(paste0("~/Documents/myeloid_subclass/training_data_CAFs/classifier_layer_4_CAF_subtypes.RData"))
top_genes_unlisted_layer_4_CAFs <- top_genes_unlisted
model_layer_4_CAFs <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_4_CAFs, overwrite = TRUE)
usethis::use_data(model_layer_4_CAFs, overwrite = TRUE)
rm(list = ls())


load(paste0("~/Documents/myeloid_subclass/classifier_outputs_myeloid/classifier_layer_5_cDC_subtypes.RData"))
top_genes_unlisted_layer_5_cDC <- top_genes_unlisted
model_layer_5_cDC <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_cDC, overwrite = TRUE)
usethis::use_data(model_layer_5_cDC, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/myeloid_subclass/classifier_outputs_myeloid/classifier_layer_5_macrophage_subtypes.RData"))
top_genes_unlisted_layer_5_macrophage <- top_genes_unlisted
model_layer_5_macrophage <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_macrophage, overwrite = TRUE)
usethis::use_data(model_layer_5_macrophage, overwrite = TRUE)
rm(list = ls())

load(paste0("~/Documents/myeloid_subclass/classifier_outputs_myeloid/classifier_layer_5_monocyte_subtypes.RData"))
top_genes_unlisted_layer_5_monocyte <- top_genes_unlisted
model_layer_5_monocyte <- rf_classifier_cell_lines
usethis::use_data(top_genes_unlisted_layer_5_monocyte, overwrite = TRUE)
usethis::use_data(model_layer_5_monocyte, overwrite = TRUE)
rm(list = ls())




