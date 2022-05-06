#' Score Class
#'
#' @param cell_name the name of the cell to score
#' @param predictions path to the matrix generated from classify cells
#' @param layer layer of the model can be layer_1, layer_2_blood, layer_2_cancer, layer_3_TNK, layer_3_MDC, layer_3_BCell, layer_4_CD4_CD8, layer_4_CD8_NK, layer_3_ov_endo, layer_3_breast_lung, layer_3_brain_lung_nbm_skin, layer_4_brain_nbm_skin, layer_5_brain_nbm
#'
#' @return A vector of the predicted class based on scores
#' @export
score_class <- function(cell_name, predictions, layer){
  layer_predictions <- predictions[cell_name,]
  if(layer == "layer_1"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "cancer_normal_stromal_score"){
      predicted_class <- "Tissue_Cell_Normal_or_Cancer"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "blood_score"){
      predicted_class <- "Blood_Cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_2_blood"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "B cell"){
      predicted_class <- "B Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "non_B_lymphocyte_score"){
      predicted_class <- "T_or_NK_lymphocyte"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "macrophage_DC_score"){
      predicted_class <- "macrophage_or_dendritic_cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "HSPC_score"){
      predicted_class <- "HSPC"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_3_TNK"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "CD4_CD8_score"){
      predicted_class <- "CD4 or CD8 T cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "NK_CD8_score"){
      predicted_class <- "NK or CD8 T cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_3_MDC"){
    if (is.na(layer_predictions)== TRUE){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "macrophage_score"){
      predicted_class <- "Macrophage"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "DC_score"){
      predicted_class <- "Dendritic Cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_3_BCell"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "B_cell_score"){
      predicted_class <- "B Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "plasmablast_score"){
      predicted_class <- "Plasmablast"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_4_CD4_CD8"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "CD4+ T cell"){
      predicted_class <- "CD4+ T cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "CD8+ T cell"){
      predicted_class <- "CD8+ T cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }

  }
  if(layer == "layer_4_CD8_NK"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "CD8+ T cell"){
      predicted_class <- "CD8+ T cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "Natural killer cell"){
      predicted_class <- "Natural killer cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_2_non_blood"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
             == "stromal_score"){
      predicted_class <- "Stromal Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "non_stromal_score"){
      predicted_class <- "Non Stromal Cell"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_3_non_stromal"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
             == "soft_tissue_neuro_score"){
      predicted_class <- "Soft Tissue or Neuro Cancer Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "non_GI_score"){
      predicted_class <- "Non GI Epithelial Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "GI_score"){
      predicted_class <- "GI Epithelial Cell"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }



  if(layer == "layer_4_non_GI"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "ov_endo_kid_score"){
      predicted_class <- "Ovarian/Endometrial/Kidney Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "breast_lung_prostate_score"){
      predicted_class <- "Breast/Lung/Prostate Cell"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_4_GI"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "billiary_score"){
      predicted_class <- "Billiary Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "digestive_score"){
      predicted_class <- "Colorectal/Esophageal/Gastric Cell"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_5_ov_endo_kid"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "Ovarian Cancer"){
      predicted_class <- "Ovarian Cancer Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "Endometrial/Uterine Cancer"){
      predicted_class <- "Endometrial Cancer Cell"
    }else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
             == "Kidney Cancer"){
      predicted_class <- "Kidney Cancer Cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_5_breast_lung_prostate"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "Breast Cancer"){
      predicted_class <- "Breast Cancer Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "Lung Cancer"){
      predicted_class <- "Lung Cancer Cell"
    }else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
             == "Prostate Cancer"){
      predicted_class <- "Prostate Cancer Cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_6_breast"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "TNBC"){
      predicted_class <- "TNBC Breast Cancer Cell"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "HER2+"){
      predicted_class <- "Her2+ Breast Cancer Cell"
    }else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
             == "ER+"){
      predicted_class <- "ER+ Breast Cancer Cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_5_biliary"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_3_stromal"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    } else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
              == "fibroblasts_score"){
      predicted_class <- "Fibroblasts"
    }else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_4_soft_tissue_neuro"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "soft_tissue_neuro_score"){
      predicted_class <- "Soft Tissue or Neuro Cancer Cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_5_soft_tissue_neuro"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
            == "brain_nbm_score"){
      predicted_class <- "Brain/Neuroblastoma Cancer Cell"
    } else if(colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
              == "soft_tissue_score"){
      predicted_class <- "Soft Tissue Cancer Cell"
    } else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_5_digestive"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_6_brain_nbm"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_6_soft_tissue"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  if(layer == "layer_4_dendritic"){
    if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
      predicted_class <- "Cell low quality"
    }
    else{
      predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions, class) == "numeric")], decreasing = T)[1])
    }
  }
  return(predicted_class)
}
