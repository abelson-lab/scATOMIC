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
                                        "Sarcoma")
    blood_classes <- c("B cell","CD4+ T cell", "CD8+ T cell", "HSPC", "Macrophage","cDC", "pDC", "ASDC","Plasmablast",
                       "Natural killer cell")
    cancer_normal_stromal_score <- rowSums(layer_predictions[,cancer_normal_stromal_classes])
    blood_score <- rowSums(layer_predictions[,blood_classes])
    layer_classes <- cbind(layer_predictions, cancer_normal_stromal_score, blood_score)
  }
  if(layer == "layer_2_blood"){
    B_cell_score <- rowSums(layer_predictions[,c("B cell", "Plasmablast")])
    non_B_lymphocyte_score <- rowSums(layer_predictions[,c("CD4+ T cell","CD8+ T cell", "Natural killer cell")])
    macrophage_DC_score <- rowSums(layer_predictions[,c("Macrophage","cDC", "pDC", "ASDC")])
    HSPC_score <- layer_predictions[,"HSPC"]
    layer_classes <- cbind(layer_predictions, B_cell_score, non_B_lymphocyte_score, macrophage_DC_score,
                           HSPC_score)
  }
  if(layer == "layer_3_TNK"){
    CD4_CD8_score <- rowSums(layer_predictions[,c("CD4+ T cell", "CD8+ T cell")])
    NK_CD8_score <- rowSums(layer_predictions[,c("Natural killer cell", "CD8+ T cell")])
    layer_classes <- cbind(layer_predictions, CD4_CD8_score, NK_CD8_score)
  }
  if(layer == "layer_3_MDC"){
    macrophage_score <- layer_predictions[,"Macrophage"]
    DC_score <- rowSums(layer_predictions[,c("ASDC", "cDC", "pDC")])
    layer_classes <- cbind(layer_predictions, macrophage_score, DC_score)
  }
  if(layer == "layer_3_BCell"){
    B_cell_score <- layer_predictions[,"B cell"]
    plasmablast_score <- layer_predictions[,"Plasmablast"]
    layer_classes <- cbind(layer_predictions, B_cell_score, plasmablast_score)
  }
  if(layer == "layer_2_non_blood"){
    stromal_score <- layer_predictions$`Endothelial Cells` +
      layer_predictions$Fibroblasts + layer_predictions$Myofibroblasts + layer_predictions$`Smooth Muscle Cells`
    non_stromal_score <- layer_predictions$`Bone Cancer`  +
      layer_predictions$Sarcoma +  layer_predictions$`Brain Cancer` +
      layer_predictions$Neuroblastoma + layer_predictions$`Skin Cancer` + layer_predictions$`Bile Duct Cancer` + layer_predictions$`Bladder Cancer` + layer_predictions$`Breast Cancer` +
      layer_predictions$`Colon/Colorectal Cancer` + layer_predictions$`Endometrial/Uterine Cancer` + layer_predictions$`Esophageal Cancer` +
      layer_predictions$`Gallbladder Cancer` + layer_predictions$`Gastric Cancer` + layer_predictions$`Kidney Cancer` + layer_predictions$`Liver Cancer` + layer_predictions$`Lung Cancer` +
      layer_predictions$`Ovarian Cancer` + layer_predictions$`Pancreatic Cancer` + layer_predictions$`Prostate Cancer`
    layer_classes <- cbind(layer_predictions,stromal_score,non_stromal_score)
  }
  if(layer == "layer_3_non_stromal"){
    non_GI_score <- layer_predictions$`Endometrial/Uterine Cancer` + layer_predictions$`Ovarian Cancer` + layer_predictions$`Breast Cancer` +
      layer_predictions$`Lung Cancer` +layer_predictions$`Prostate Cancer` + layer_predictions$`Kidney Cancer`
    GI_score <- layer_predictions$`Bile Duct Cancer` + layer_predictions$`Bladder Cancer` + layer_predictions$`Colon/Colorectal Cancer` +
      layer_predictions$`Esophageal Cancer` + layer_predictions$`Gallbladder Cancer` + layer_predictions$`Gastric Cancer` +
      layer_predictions$`Liver Cancer` + layer_predictions$`Pancreatic Cancer`
    soft_tissue_neuro_score <- layer_predictions$`Bone Cancer` +
      layer_predictions$Sarcoma +  layer_predictions$`Brain Cancer` +
      layer_predictions$Neuroblastoma + layer_predictions$`Skin Cancer` + layer_predictions$`Lung Cancer`
    layer_classes <- cbind(layer_predictions,soft_tissue_neuro_score, non_GI_score, GI_score)
  }
  if(layer == "layer_4_non_GI"){
    ov_endo_kid_score <- layer_predictions$`Endometrial/Uterine Cancer` + layer_predictions$`Ovarian Cancer`  + layer_predictions$`Kidney Cancer`
    breast_lung_prostate_score <- layer_predictions$`Breast Cancer` + layer_predictions$`Lung Cancer` +layer_predictions$`Prostate Cancer`
    layer_classes <- cbind(layer_predictions, ov_endo_kid_score, breast_lung_prostate_score)
  }
  if(layer == "layer_4_GI"){
    digestive_score <- layer_predictions$`Colon/Colorectal Cancer` + layer_predictions$`Esophageal Cancer` + layer_predictions$`Gastric Cancer`
    billiary_score <- layer_predictions$`Bile Duct Cancer` + layer_predictions$`Bladder Cancer` + layer_predictions$`Gallbladder Cancer` +
      layer_predictions$`Liver Cancer` + layer_predictions$`Pancreatic Cancer`
    layer_classes <- cbind(layer_predictions, billiary_score, digestive_score)
  }
  if(layer == "layer_4_soft_tissue_neuro"){
    soft_tissue_neuro_score <- layer_predictions$`Brain Cancer` + layer_predictions$Neuroblastoma +
      layer_predictions$`Bone Cancer`  + layer_predictions$Sarcoma
    layer_classes <- cbind(layer_predictions, soft_tissue_neuro_score)
  }
  if(layer == "layer_5_soft_tissue_neuro"){
    brain_nbm_score <- layer_predictions$`Brain Cancer` + layer_predictions$Neuroblastoma
    soft_tissue_score <-   layer_predictions$`Bone Cancer`  + layer_predictions$Sarcoma
    layer_classes <- cbind(layer_predictions, brain_nbm_score, soft_tissue_score)
  }
  if(layer == "layer_3_stromal"){
    fibroblasts_score <-   layer_predictions$Fibroblasts + layer_predictions$Myofibroblasts
    layer_classes <- cbind(layer_predictions, fibroblasts_score)
  }
  if(layer %in% c("layer_4_dendritic","layer_5_biliary",
                  "layer_6_soft_tissue", "layer_6_brain_nbm",
                  "layer_5_digestive", "layer_5_breast_lung_prostate",
                  "layer_5_ov_endo_kid", "layer_4_CD4_CD8",
                  "layer_4_CD8_NK", "layer_6_breast")){
    print("nothing to score in this layer")
    layer_classes <- layer_predictions
  }
  return(layer_classes)
}
