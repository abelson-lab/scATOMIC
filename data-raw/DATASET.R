## code to prepare `DATASET` dataset goes here
load("~/Documents/classifier_tests_nov21/classifiers_final_0.4/nov22_blood_corrected_layer_1_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_1 <- top_genes_unlisted
model_layer_1 <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_1, overwrite = TRUE)
usethis::use_data(model_layer_1, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifier_tests_nov21/classifiers_final_0.4/nov22_blood_corrected_layer_2_blood_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_2_blood <- top_genes_unlisted
model_layer_2_blood <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_2_blood, overwrite = TRUE)
usethis::use_data(model_layer_2_blood, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_2_non_blood_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_2_normal_tissue_cancer <- top_genes_unlisted
model_layer_2_normal_tissue_cancer <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_2_normal_tissue_cancer, overwrite = TRUE)
usethis::use_data(model_layer_2_normal_tissue_cancer, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_3_B_cell_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_3_BCell <- top_genes_unlisted
model_layer_3_BCell <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_BCell, overwrite = TRUE)
usethis::use_data(model_layer_3_BCell, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifier_tests_nov21/classifiers_final_0.4/nov22_blood_corrected_layer_3_MDC_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_3_MDC <- top_genes_unlisted
model_layer_3_MDC <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_MDC, overwrite = TRUE)
usethis::use_data(model_layer_3_MDC, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_3_non_stromal_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_3_non_stromal <- top_genes_unlisted
model_layer_3_non_stromal <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_non_stromal, overwrite = TRUE)
usethis::use_data(model_layer_3_non_stromal, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_3_stromal_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_3_stromal <- top_genes_unlisted
model_layer_3_stromal <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_stromal, overwrite = TRUE)
usethis::use_data(model_layer_3_stromal, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_3_TNK_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_3_TNK <- top_genes_unlisted
model_layer_3_TNK <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_3_TNK, overwrite = TRUE)
usethis::use_data(model_layer_3_TNK, overwrite = TRUE)
rm(list = ls())


load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_4_CD4_CD8_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_4_CD4_CD8 <- top_genes_unlisted
model_layer_4_CD4_CD8 <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_CD4_CD8, overwrite = TRUE)
usethis::use_data(model_layer_4_CD4_CD8, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_4_CD8_NK_cell_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_4_CD8_NK <- top_genes_unlisted
model_layer_4_CD8_NK <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_CD8_NK, overwrite = TRUE)
usethis::use_data(model_layer_4_CD8_NK, overwrite = TRUE)
rm(list = ls())


load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_4_DC_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_4_DC_cell <- top_genes_unlisted
model_layer_4_DC_cell <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_DC_cell, overwrite = TRUE)
usethis::use_data(model_layer_4_DC_cell, overwrite = TRUE)
rm(list = ls())


load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_4_GI_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_4_GI <- top_genes_unlisted
model_layer_4_GI <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_GI, overwrite = TRUE)
usethis::use_data(model_layer_4_GI, overwrite = TRUE)
rm(list = ls())


load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_4_non_GI_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_4_non_GI <- top_genes_unlisted
model_layer_4_non_GI <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_non_GI, overwrite = TRUE)
usethis::use_data(model_layer_4_non_GI, overwrite = TRUE)
rm(list = ls())


load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_4_soft_tissue_neuro_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_4_soft_neuro_cancer <- top_genes_unlisted
model_layer_4_soft_neuro_cancer <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_4_soft_neuro_cancer, overwrite = TRUE)
usethis::use_data(model_layer_4_soft_neuro_cancer, overwrite = TRUE)
rm(list = ls())


load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_5_biliary_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_5_biliary <- top_genes_unlisted
model_layer_5_biliary <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_biliary, overwrite = TRUE)
usethis::use_data(model_layer_5_biliary, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_5_breast_lung_prostate_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_5_breast_lung_prostate <- top_genes_unlisted
model_layer_5_breast_lung_prostate <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_breast_lung_prostate, overwrite = TRUE)
usethis::use_data(model_layer_5_breast_lung_prostate, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_5_digestive_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_5_digestive <- top_genes_unlisted
model_layer_5_digestive <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_digestive, overwrite = TRUE)
usethis::use_data(model_layer_5_digestive, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_5_ov_endo_kid_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_5_ov_endo_kid <- top_genes_unlisted
model_layer_5_ov_endo_kid  <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_ov_endo_kid, overwrite = TRUE)
usethis::use_data(model_layer_5_ov_endo_kid, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_5_soft_tissue_neuro_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin <- top_genes_unlisted
model_layer_5_soft_neuro_cancer_no_lung_skin <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin, overwrite = TRUE)
usethis::use_data(model_layer_5_soft_neuro_cancer_no_lung_skin, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_6_brain_nbm_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_6_brain_nbm <- top_genes_unlisted
model_layer_6_brain_nbm <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_6_brain_nbm, overwrite = TRUE)
usethis::use_data(model_layer_6_brain_nbm, overwrite = TRUE)
rm(list = ls())

load("~/Documents/classifiers_testing_paper_july2021/classifiers_final_0.4/oct7_layer_6_soft_tissue_0_sds_min50_max200_pct2_0.4.RData")
top_genes_unlisted_layer_6_soft_tissue <- top_genes_unlisted
model_layer_6_soft_tissue <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_6_soft_tissue, overwrite = TRUE)
usethis::use_data(model_layer_6_soft_tissue, overwrite = TRUE)
rm(list = ls())


load("~/Documents/classifiers_testing_paper_july2021/breast_subclassifier/oct19_layer_6_breast_subclassifier_0_sds_0.4_pct.2.RData")
top_genes_unlisted_layer_6_breast <- top_genes_unlisted
model_layer_6_breast <- rf_classifier_cell_lines

usethis::use_data(top_genes_unlisted_layer_6_breast, overwrite = TRUE)
usethis::use_data(model_layer_6_breast, overwrite = TRUE)
rm(list = ls())

demo_lung_data <- readRDS("~/Documents/lung_demo_data.RDS")
usethis::use_data(demo_lung_data, overwrite = TRUE)
rm(list = ls())

list_markers <- list.files("~/Downloads/expression_full/", full.names = T)
markers_master <- data.frame()
for(i in 1:length(list_markers)){
  markers_file <- fread(list_markers[i])
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












