#' Add Unclassified
#'
#' @param cell_name the name of the cell to score
#' @param predictions path to the matrix generated from classify cells and merged with score class
#' @param layer layer of the model can be layer_1, layer_2_blood, layer_2_cancer, layer_3_TNK, layer_3_MDC, layer_3_BCell, layer_4_CD4_CD8, layer_4_CD8_NK, layer_3_ov_endo, layer_3_breast_lung, layer_3_brain_lung_nbm_skin, layer_4_brain_nbm_skin, layer_5_brain_nbm
#' @param cutoff value for cutoff
#'
#' @return vecotr of filtered classes
#' @export
add_unclassified_automatic_threshold <- function(cell_name, predictions, layer, threshold_use){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  layer_predictions <- predictions[cell_name,]

  if(layer == "layer_1"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- layer_predictions["predicted_class"]
      } else {
        predicted_with_cutoff <- "unclassified_any_cell"
      }
    }
  }
  if(layer == "layer_2_non_blood"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- layer_predictions["predicted_class"]
      } else {
        predicted_with_cutoff <- "unclassified_normal_or_cancer_tissue"
      }
    }}
  if(layer == "layer_2_blood"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- layer_predictions["predicted_class"]
      }
      else {
        predicted_with_cutoff <- "unclassified_blood_cell"
      }
    }}
  if(layer == "layer_3_TNK"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- layer_predictions["predicted_class"]
      }
      else {
        predicted_with_cutoff <- "unclassified_T_or_NK_cell"
      }
    }}
  if(layer == "layer_3_MDC"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- layer_predictions[ "predicted_class"]
      }
      else {
        predicted_with_cutoff <- "unclassified_macrophage_or_DC"
      }
    }}
  if(layer == "layer_3_BCell"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "unclassified_B_cell_or_plasmablast"
      }
    }}
  if(layer == "layer_4_CD4_CD8"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- layer_predictions[ "predicted_class"]
      }
      else {
        predicted_with_cutoff <- "CD4 or CD8 T cell"
      }
    }}
  if(layer == "layer_4_CD8_NK"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- layer_predictions[ "predicted_class"]
      }
      else {
        predicted_with_cutoff <- "CD8 T or NK cell"
      }
    }}
  if(layer == "layer_4_non_GI"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "unclassified_non_GI_epithelial_cell"
      }
    }}
  if(layer == "layer_4_GI"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "unclassified_GI_epithelial_cell"
      }
    }}
  if(layer == "layer_3_stromal"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Stromal Cell"
      }
    }}
  if(layer == "layer_3_non_stromal"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Non Stromal Cell"
      }
    }}
  if(layer == "layer_4_soft_tissue_neuro"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Unclassified Soft Tissue or Neuro Cancer Cell"
      }
    }}
  if(layer == "layer_5_soft_tissue_neuro"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Soft Tissue or Neuro Cancer Cell"
      }
    }}

  if(layer == "layer_5_ov_endo_kid"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Ovarian/Endometrial/Kidney"
      }
    }}
  if(layer == "layer_5_breast_lung_prostate"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Breast/Lung/Prostate"
      }
    }}
  if(layer == "layer_6_breast"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Breast Cancer Cell"
      }
    }}
  if(layer == "layer_5_biliary"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Biliary/Hepatic Cancer Cell"
      }
    }}
  if(layer == "layer_5_digestive"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Digestive Tract Cancer Cell"
      }
    }}
  if(layer == "layer_6_brain_nbm"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Brain/Neuroblastoma Cancer Cell"
      }
    }}
  if(layer == "layer_6_soft_tissue"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Soft Tissue Cancer Cell"
      }
    }}
  if(layer == "layer_4_dendritic"){
    if ((layer_predictions["predicted_class"]) == "Cell low quality"){
      predicted_with_cutoff <- "Cell low quality"
    } else {
      scores <- which(lapply(layer_predictions, class) == "numeric")
      max_score <- max(layer_predictions[which(lapply(layer_predictions, class) == "numeric")])
      max_score_class <- colnames(layer_predictions)[which(layer_predictions == max_score)]
      if(length(max_score_class) > 1){
        max_score_class <- max_score_class[length(max_score_class)]
      }
      if(max_score > threshold_use[[max_score_class]]){
        predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
      }
      else {
        predicted_with_cutoff <- "Dendritic Cell"
      }
    }}
  return(predicted_with_cutoff)
}
