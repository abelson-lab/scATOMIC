#' outersect
#'
#' @param x first element
#' @param y second element
#'
#' @return vector of non overlapping features - opposite of interesect
#' @export
outersect <- function(x, y) {
  non_overlapping <- sort(c(setdiff(x, y),
                            setdiff(y, x)))
  return(non_overlapping)
}
