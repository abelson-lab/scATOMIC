% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_new_scATOMIC_layer.R
\name{classify_new_scATOMIC_layer}
\alias{classify_new_scATOMIC_layer}
\title{classify_new_scATOMIC_layer}
\usage{
classify_new_scATOMIC_layer(
  rf_model,
  selected_features,
  rna_counts,
  cell_names,
  layer_name,
  mc.cores
)
}
\arguments{
\item{rf_model}{model trained by get_new_scATOMIC_layer()}

\item{selected_features}{selected features returned by get_new_scATOMIC_layer()}

\item{rna_counts}{gene by cell count matrix of query dataset}

\item{cell_names}{names of cells to be classified in layer. We recommend passing row.names of the output of create_summary_matrix() of cells that are to be further classified}

\item{layer_name}{name of new layer, unclassified cells will receive the name Unclassified_Cell_from_layer_name}

\item{mc.cores}{number of computing cores to use for parallel processing}
}
\value{
returns a matrix containing cell classifications and each class score
}
\description{
classify_new_scATOMIC_layer
}
\examples{
\dontrun{
#visit github for more information on how to generate training datasets.
breast_cancer_subclassification <- get_new_scATOMIC_layer(training_data = Wu_et_al_breast_count_mat,cell_type_metadata = Wu_et_al_breast_metadata$subtype,
output_dir = "Wu_etal_2021_BRCA_scRNASeq/breast_subclassification_layer/",
mc.cores = 6, layer_name = "breast_subclassification", n_cells_replicate = 10000, n_trees = 500)
predictions_Pal_0125 <- run_scATOMIC(Pal_0125)
results_Pal_0125 <- create_summary_matrix(prediction_list = predictions_Pal_0125)
cells_to_subclassify <- row.names(results_Pal_0125)[which(results_Pal_0125$scATOMIC_pred == "Breast Cancer Cell")]
breast_subclassifications <- classify_new_scATOMIC_layer(rf_model = breast_cancer_subclassification[["rf_classifier"]], selected_features = breast_cancer_subclassification[["selected_features"]],
                                                   rna_counts = Pal_0125, cell_names = cells_to_subclassify, layer_name = "Breast Cancer Cells", mc_cores = 6  )


}

}
