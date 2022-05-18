#' Class for Cutoff
#'
#' @param cell_name the name of the cell
#' @param predictions the dataframe of predictions
#'
#' @return class to use in automatic threshold function
#' @export
#'
class_for_cutoff <- function(cell_name, predictions){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  layer_predictions <- predictions[cell_name,]
  if(levels(as.factor(is.na(layer_predictions) %in% TRUE))){
    class_to_use <- NA
  } else{
    scores <- which(lapply(layer_predictions, class) == "numeric")
    class_to_use <- colnames(layer_predictions)[which(layer_predictions == max(layer_predictions[scores]))]
    if(length(class_to_use) > 1){
      class_to_use <- class_to_use[length(class_to_use)]
    }
  }
  return(class_to_use)
}
