#' classify_layer
#'
#' @param rna_counts count matrix
#' @param cells_to_use vector indicating which cells to classify in layer
#' @param imputation whether to apply MAGIC imputation - recommended
#' @param genes_in_model genes used in the model
#' @param model model
#' @param ref_based whether to use ref based - currently depreceated
#' @param mc.cores number of cores
#' @param unimodal_nsd number of sd for setting threshold in unimodal distributions
#' @param bimodal_nsd number of sd for setting threshold in bimodal distributions
#' @param layer which layer is being classified
#'
#' @return returns a prediction matrix for a layer - used in run_scATOMIC
#' @export

classify_layer <- function(rna_counts, cells_to_use, imputation = TRUE, genes_in_model, model,
                           ref_based = F, mc.cores = (parallel::detectCores()-1), unimodal_nsd = 3, bimodal_nsd = 2,
                           layer, normalized_counts){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  layer_predictions <- classify_cells_RNA_no_scale(rna_counts = rna_counts, imputation = imputation,
                                                                       genes_in_model = genes_in_model,
                                                                       model = model, cells_to_use = cells_to_use,
                                                                       ref_based = ref_based, normalized_counts=normalized_counts, mc.cores = mc.cores)
  layer_predictions <- add_class_per_layer(layer_predictions = layer_predictions,
                                                               layer = layer)

  predicted_class <- unlist(parallel::mclapply(row.names(layer_predictions),
                                               score_class, predictions = layer_predictions,
                                               layer = layer, mc.cores = mc.cores), use.names = F)
  layer_predictions <- cbind(layer_predictions, predicted_class)
  layer_predictions$predicted_class <- as.character(predicted_class)
  threshold_per_class <- get_auto_threshold(layer_predictions, mc.cores=mc.cores)
  predicted_tissue_with_cutoff <- unlist(parallel::mclapply(row.names(layer_predictions),
                                                            add_unclassified_automatic_threshold, predictions = layer_predictions,
                                                            layer = layer, threshold_use = threshold_per_class,
                                                            mc.cores = mc.cores), use.names = F)
  layer_predictions <- cbind(layer_predictions, predicted_tissue_with_cutoff)
  return(layer_predictions)
}
