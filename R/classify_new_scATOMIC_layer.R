#' classify_new_scATOMIC_layer
#'
#' @param rf_model model trained by get_new_scATOMIC_layer()
#' @param selected_features selected features returned by get_new_scATOMIC_layer()
#' @param rna_counts gene by cell count matrix of query dataset
#' @param cell_names names of cells to be classified in layer. We recommend passing row.names of the output of create_summary_matrix() of cells that are to be further classified
#' @param layer_name name of new layer, unclassified cells will receive the name Unclassified_Cell_from_layer_name
#' @param mc_cores number of computing cores to use for parallel processing
#'
#' @return returns a matrix containing cell classifications and each class score
#' @export
#' @examples
#' \dontrun{
#' #visit github for more information on how to generate training datasets.
#' breast_cancer_subclassification <- get_new_scATOMIC_layer(training_data = Wu_et_al_breast_count_mat,cell_type_metadata = Wu_et_al_breast_metadata$subtype,
#' output_dir = "Wu_etal_2021_BRCA_scRNASeq/breast_subclassification_layer/",
#' mc_cores = 6, layer_name = "breast_subclassification", n_cells_replicate = 10000, n_trees = 500)
#' predictions_Pal_0125 <- run_scATOMIC(Pal_0125)
#' results_Pal_0125 <- create_summary_matrix(prediction_list = predictions_Pal_0125)
#' cells_to_subclassify <- row.names(results_Pal_0125)[which(results_Pal_0125$scATOMIC_pred == "Breast Cancer Cell")]
#' breast_subclassifications <- classify_new_scATOMIC_layer(rf_model = breast_cancer_subclassification[["rf_classifier"]], selected_features = breast_cancer_subclassification[["selected_features"]],
#'                                                    rna_counts = Pal_0125, cell_names = cells_to_subclassify, layer_name = "Breast Cancer Cells", mc_cores = 6  )
#'
#'
#' }
#'
classify_new_scATOMIC_layer <- function(rf_model, selected_features, rna_counts, cell_names, layer_name, mc_cores){
  normalized_counts <- Rmagic::library.size.normalize(t(as.matrix(rna_counts)))
  normalized_counts <- t(sqrt(normalized_counts))
  layer_predictions <- scATOMIC::classify_cells_RNA_no_scale(rna_counts =rna_counts , imputation = T,
                                                             genes_in_model = selected_features,
                                                             model = rf_model, cells_to_use = cell_names,
                                                             ref_based = F, normalized_counts=normalized_counts)
  predicted_class <- unlist(parallel::mclapply(row.names(layer_predictions),
                                               function (cell_name, predictions)
                                               {    layer_predictions <- predictions[cell_name, ]
                                               predicted_class <- colnames(sort(layer_predictions[which(lapply(layer_predictions,
                                                                                                               class) == "numeric")], decreasing = T)[1])}, predictions = layer_predictions,
                                               mc.cores = mc_cores), use.names = F)
  layer_predictions <- cbind(layer_predictions, predicted_class)
  layer_predictions$predicted_class <- as.character(predicted_class)
  threshold_per_class <- scATOMIC::get_auto_threshold(layer_predictions)
  predicted_tissue_with_cutoff <- unlist(parallel::mclapply(row.names(layer_predictions),
                                                            function (cell_name, predictions, threshold_use, layer_name)
                                                            {
                                                              layer_predictions <- predictions[cell_name, ]
                                                              scores <- which(lapply(layer_predictions, class) ==
                                                                                "numeric")
                                                              max_score <- max(layer_predictions[which(lapply(layer_predictions,
                                                                                                              class) == "numeric")])
                                                              max_score_class <- colnames(layer_predictions)[which(layer_predictions ==
                                                                                                                     max_score)]
                                                              if (length(max_score_class) > 1) {
                                                                max_score_class <- max_score_class[length(max_score_class)]
                                                              }
                                                              if (max_score > threshold_use[[max_score_class]]) {
                                                                predicted_with_cutoff <- as.character(layer_predictions["predicted_class"])
                                                              }
                                                              else {
                                                                predicted_with_cutoff <- paste0("Unclassified_Cell_from_", layer_name)
                                                              }
                                                            }, predictions = layer_predictions,
                                                            threshold_use = threshold_per_class, layer_name = layer_name,
                                                            mc.cores = 6), use.names = F)
  layer_predictions <- cbind(layer_predictions, predicted_tissue_with_cutoff)
  return(layer_predictions)
}




