#' get_auto_threshold
#'
#' @param predictions dataframe of scores and predicted class
#' @param mc.cores number of cores for parallel computing
#' @param unimodal_nsd number of standard deviations if scores are unimodal
#' @param bimodal_nsd number of standard deviations if scores are multimodal
#'
#' @return returns a confidence threshold for each cell type class
#' @export
get_auto_threshold <- function(predictions, mc.cores = (parallel::detectCores()-1), unimodal_nsd = 3, bimodal_nsd = 2){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  cutoff_class <- unlist(parallel::mclapply(row.names(predictions), scATOMIC::class_for_cutoff, predictions = predictions, mc.cores = mc.cores), use.names = F)
  predictions <- cbind(na.omit(predictions), na.omit(cutoff_class))
  scores_to_get_threshold <- levels(as.factor(cutoff_class))
  list_scores_for_cutoffs <- list()
  for(i in 1:length(scores_to_get_threshold)){
    list_scores_for_cutoffs[[i]] <- predictions[,scores_to_get_threshold[i]]
    names(list_scores_for_cutoffs)[i] <- scores_to_get_threshold[i]
  }
  threshold_per_class <- lapply(list_scores_for_cutoffs, scATOMIC::automatic_threshold, unimodal_nsd = unimodal_nsd, bimodal_nsd = bimodal_nsd)
  return(threshold_per_class)
}
