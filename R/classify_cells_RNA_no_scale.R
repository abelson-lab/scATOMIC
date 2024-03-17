#' classify_cells_RNA_no_scale
#'
#' @param imputation TRUE or FALSE for whether to perform imputation
#' @param model Path to RF model to use
#' @param genes_in_model Vector of genes to use in the model
#' @param cells_to_use Vector of cell names to use from the rna_counts, refers to their column names
#' @param normalized_counts library-size normalized counts
#' @param mc.cores number of cores
#'
#' @return A matrix with predictions of each cell class in the model
#' @export
classify_cells_RNA_no_scale <- function( imputation, model, genes_in_model, cells_to_use, normalized_counts, mc.cores = (parallel::detectCores()-1)){
  #normalize via magic
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  if(length(cells_to_use)==0){
    print("There are no cells to be classified in this layer")
  } else {

    #get genes shared between the model and input database
    shared_genes <- dplyr::intersect(row.names(normalized_counts), genes_in_model)
    if(length(cells_to_use) == 1){
      filtered_counts <- as.data.frame(na.omit(normalized_counts[shared_genes,cells_to_use]))
      colnames(filtered_counts) <- cells_to_use
    } else {
      filtered_counts <- na.omit(normalized_counts[shared_genes,cells_to_use])
    }
    NaNs <- which(filtered_counts[, 1] == "NaN")
    if(length(NaNs > 0)){
      filtered_counts <- filtered_counts[-NaNs, ]
    }
    #run magic imputation or not
    if(imputation == TRUE & length(cells_to_use) >= 20){
      filtered_counts <- t(as.matrix(filtered_counts))
      defaultW <- getOption("warn")
      options(warn = -1)
      sc_magic <- Rmagic::magic(filtered_counts,verbose = F, seed=123)
      options(warn = defaultW)
      filtered_counts <- as.data.frame(t(as.data.frame(sc_magic)))

    }
    if(imputation == FALSE | length(cells_to_use) < 10){
      filtered_counts <- as.data.frame(filtered_counts)
    }

    outersect <- function(x, y) {
      sort(c(setdiff(x, y),
             setdiff(y, x)))
    }
    missing_genes <- outersect(genes_in_model, shared_genes)
    if(length(NaNs > 0)){
      missing_genes <- c(missing_genes, names(NaNs))
    }
    #fill in missing features as 0
    if(length(missing_genes) > 0){
      missing_genes_mat <- as.data.frame(matrix(0, nrow = length(missing_genes),
                                                ncol = ncol(filtered_counts)))
      row.names(missing_genes_mat) <- missing_genes
      colnames(missing_genes_mat) <- colnames(filtered_counts)
      filtered_counts <- rbind(filtered_counts, missing_genes_mat)
    }
    #get scaled, percent of reads per cell
    counts_as_fraction <- as.data.frame(apply(filtered_counts, 2, function(x){
      x/sum(x)
    }))
    scaled_counts_as_fraction <- as.data.frame(t(counts_as_fraction))
    colnames(scaled_counts_as_fraction) <- gsub("-", "_", colnames(scaled_counts_as_fraction))
    #apply model
    predictions_probability <- as.data.frame(predict(model, scaled_counts_as_fraction, type = "prob"))
    return(predictions_probability)

  }
}
