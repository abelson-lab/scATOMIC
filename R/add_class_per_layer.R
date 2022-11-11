#' Add Class Per Layer
#'
#' @param layer_predictions The matrix of predictions from classify_cells
#' @param layer layer of the model can be layer_1, layer_2_blood, layer_2_cancer, layer_3_TNK, layer_3_MDC, layer_3_BCell, layer_4_CD4_CD8, layer_4_CD8_NK, layer_3_ov_endo, layer_3_breast_lung, layer_3_brain_lung_nbm_skin, layer_4_brain_nbm_skin, layer_5_brain_nbm
#'
#' @return matrix with classes
#' @export
#'
add_class_per_layer <- function(layer_predictions, layer){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  if(layer == "layer_1"){
    cancer_normal_stromal_classes <- c( "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer",     "Brain Cancer",    "Breast Cancer",
                                        "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
                                        "Gallbladder Cancer", "Gastric Cancer", "Glial Cells",  "Kidney Cancer", "Liver Cancer", "Lung Cancer",
                                        "Neuroblastoma",  "Oligodendrocytes","Ovarian Cancer",
                                        "Pancreatic Cancer",
                                        "Prostate Cancer", "Skin Cancer", "Endothelial Cells", "Fibroblasts", "Myofibroblasts",  "Smooth Muscle Cells",
                                        "Sarcoma","Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts")
    blood_classes <- c("ASDC","B cell","CD4+ T cell", "CD8+ T cell", "Macrophage","Monocyte","Mast cell","cDC", "pDC","Plasmablast",
                       "Natural killer cell", "HSPC")
    cancer_normal_stromal_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% cancer_normal_stromal_classes)]))
    blood_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% blood_classes)]))
    layer_classes <- cbind(layer_predictions, cancer_normal_stromal_score, blood_score)
  }
  if(layer == "layer_2_blood"){
    B_cell_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("B cell", "Plasmablast"))]))
    non_B_lymphocyte_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in%c("CD4+ T cell","CD8+ T cell", "Natural killer cell"))]))
    macrophage_DC_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in%c("ASDC","Macrophage","cDC", "pDC", "Monocyte"))]))
    mast_score <- layer_predictions[,which(colnames(layer_predictions) %in% "Mast cell")]
    HSPC_score <- layer_predictions[,which(colnames(layer_predictions) %in% "HSPC")]
    layer_classes <- cbind(layer_predictions, B_cell_score, non_B_lymphocyte_score, macrophage_DC_score,
                           mast_score,HSPC_score)
  }
  if(layer == "layer_3_TNK"){
    CD4_CD8_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in%c("CD4+ T cell", "CD8+ T cell"))]))
    NK_CD8_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in%c("Natural killer cell", "CD8+ T cell"))]))
    layer_classes <- cbind(layer_predictions, CD4_CD8_score, NK_CD8_score)
  }
  if(layer == "layer_3_MDC"){
    macrophage_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Macrophage", "Monocyte"))]))
    DC_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in%c( "ASDC","cDC", "pDC"))]))
    layer_classes <- cbind(layer_predictions, macrophage_score, DC_score)
  }
  if(layer == "layer_3_BCell"){
    B_cell_score <- layer_predictions[,which(colnames(layer_predictions) == "B cell")]
    plasmablast_score <- layer_predictions[,which(colnames(layer_predictions) == "Plasmablast")]
    layer_classes <- cbind(layer_predictions, B_cell_score, plasmablast_score)
  }
  if(layer == "layer_2_non_blood"){
    stromal_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Endothelial Cells",
                                                                                                       "Fibroblasts", "Myofibroblasts", "Smooth Muscle Cells",
                                                                                                       "Cancer Associated Fibroblasts", "Cancer Associated Myofibroblasts"))]))
    non_stromal_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Bone Cancer",
                                                                                                           "Sarcoma","Brain Cancer",
                                                                                                           "Neuroblastoma", "Skin Cancer","Bile Duct Cancer", "Bladder Cancer", "Breast Cancer",
                                                                                                           "Colon/Colorectal Cancer", "Endometrial/Uterine Cancer", "Esophageal Cancer",
                                                                                                           "Gallbladder Cancer", "Gastric Cancer",  "Kidney Cancer", "Liver Cancer",  "Lung Cancer",
                                                                                                           "Ovarian Cancer", "Pancreatic Cancer", "Prostate Cancer"))]))
    layer_classes <- cbind(layer_predictions,stromal_score,non_stromal_score)
  }
  if(layer == "layer_3_non_stromal"){
    non_GI_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Breast Cancer",
                                                                                                      "Endometrial/Uterine Cancer", "Kidney Cancer",   "Lung Cancer",
                                                                                                      "Ovarian Cancer", "Prostate Cancer"))]))

    GI_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Bile Duct Cancer", "Bladder Cancer",
                                                                                                  "Colon/Colorectal Cancer", "Esophageal Cancer",
                                                                                                  "Gallbladder Cancer", "Gastric Cancer", "Liver Cancer",
                                                                                                  "Pancreatic Cancer"))]))

    soft_tissue_neuro_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Bone Cancer",
                                                                                                                 "Sarcoma","Brain Cancer",
                                                                                                                 "Neuroblastoma", "Skin Cancer", "Lung Cancer"))]))
    layer_classes <- cbind(layer_predictions,soft_tissue_neuro_score, non_GI_score, GI_score)
  }
  if(layer == "layer_4_non_GI"){
    ov_endo_kid_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Endometrial/Uterine Cancer", "Kidney Cancer",
                                                                                                           "Ovarian Cancer"))]))
    breast_lung_prostate_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Breast Cancer",
                                                                                                                    "Lung Cancer", "Prostate Cancer"))]))

    layer_classes <- cbind(layer_predictions, ov_endo_kid_score, breast_lung_prostate_score)
  }
  if(layer == "layer_4_GI"){
    digestive_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c(
      "Colon/Colorectal Cancer", "Esophageal Cancer",
      "Gastric Cancer"))]))

    billiary_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Bile Duct Cancer", "Bladder Cancer",
                                                                                                        "Gallbladder Cancer", "Liver Cancer",
                                                                                                        "Pancreatic Cancer"))]))

    layer_classes <- cbind(layer_predictions, billiary_score, digestive_score)
  }
  if(layer == "layer_4_soft_tissue_neuro"){
    soft_tissue_neuro_score <-rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Bone Cancer",
                                                                                                                "Sarcoma","Brain Cancer",
                                                                                                                "Neuroblastoma"))]))
    layer_classes <- cbind(layer_predictions, soft_tissue_neuro_score)
  }
  if(layer == "layer_5_soft_tissue_neuro"){
    brain_nbm_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Brain Cancer",
                                                                                                         "Neuroblastoma"))]))
    soft_tissue_score <-   rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Bone Cancer",
                                                                                                             "Sarcoma"))]))
    layer_classes <- cbind(layer_predictions, brain_nbm_score, soft_tissue_score)
  }
  if(layer == "layer_3_stromal"){
    fibroblasts_score <-   rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c(
      "Fibroblasts", "Myofibroblasts"))]))
    cancer_associated_fibroblast_score <- rowSums(as.data.frame(layer_predictions[,which(colnames(layer_predictions) %in% c("Cancer Associated Fibroblasts", "Cancer Associated Myofibroblasts"))]))

    layer_classes <- cbind(layer_predictions, fibroblasts_score,cancer_associated_fibroblast_score)
  }
  if(layer %in% c("layer_4_macrophage","layer_4_dendritic","layer_5_biliary",
                  "layer_6_soft_tissue", "layer_6_brain_nbm",
                  "layer_5_digestive", "layer_5_breast_lung_prostate",
                  "layer_5_ov_endo_kid", "layer_4_CD4_CD8","layer_4_CD8_NK",
                  "layer_6_breast", "layer_5_CD4", "layer_5_CD8")){
    print("nothing to score in this layer")
    layer_classes <- layer_predictions
  }
  return(layer_classes)
}
