#' Automate Threshold
#'
#' @param score Vector containing scores of a class
#' @param unimodal_nsd number of SDs to take from mean if unimodal/one-class
#' @param bimodal_nsd number of SDs to take from mean if bimodal/multi-class
#'
#' @return Value to use as a cutoff
#' @export
automatic_threshold <- function(score, unimodal_nsd, bimodal_nsd){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  score <- na.omit(score)
  if(length(score) < 50 | min(score) > 0.8){
    threshold = (min(score) - 0.01)
  } else{
    rounded_score <- round(score, digits = 1)
    rounded_collapsed_freq_vector <- agrmt::collapse(rounded_score)
    if(length(rounded_collapsed_freq_vector) < 3){
      rounded_collapsed_freq_vector <- agrmt::collapse(score)
    }
    modality <- agrmt::ajus(rounded_collapsed_freq_vector)[["type"]]
    if (modality == "A"|modality == "J"|modality == "L"|modality == "F"){
      threshold <- mean(score) - unimodal_nsd*sd(score)
    } else if(modality == "S"|modality == "U"){
      score <- as.numeric(gsub("^0$", 0.01, score))
      #index_1 <- which(score == 1)
      index_1 <- which(score > 0.89)
      index_2 <- which(score < 0.11)
      if(length(index_1) > 0){
        #score[index_1[1:(round(length(index_1))/2)]] <- as.numeric(gsub("^1$", 0.98, score[index_1[1:(round(length(index_1))/2)]]))
        #score[index_1[((round(length(index_1))/2)+1):length(index_1)]] <- as.numeric(gsub("^1$", 0.99, score[index_1[((round(length(index_1))/2)+1):length(index_1)]]))
        score[index_1] <- as.numeric(runif(length(index_1), min = 0.9, max = 1))

      }
      if(length(index_2) > 0){
        score[index_2] <- as.numeric(runif(length(index_2), min = 0.01, max = 0.1))

      }

      mixed_model <- cutoff::em(score,"normal","normal")
      mu1 <- mixed_model$param[1]
      sigma1 <- mixed_model$param[2]
      mu2 <- mixed_model$param[3]
      sigma2 <- mixed_model$param[4]
      if(mu1 > 0.5){
        threshold <- mu1 - bimodal_nsd*sigma1
      } else{
        threshold <- mu2 - bimodal_nsd*sigma2
      }
      names(threshold) <- NULL
    }
    if(threshold > 0.7){
      threshold <- 0.7
    }
    if(threshold < 0){
      threshold <- 0.6
    }
  }

  return(threshold)
}
