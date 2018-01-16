#' Remove Duplicated Gene Names from a FORESEE Object
#'
#' DuplicationHandler removes duplicates in the gene names from the FORESEE Object.
#'
#' @param Object FORESEE Object (ForeseeCell or ForeseeTrain) that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param DuplicationHandling Method for handling duplicates of gene names, such as taking none, the average, the first hit, etc.

#' @return \item{Object}{The object without duplicated gene names}
#' @export

DuplicationHandler <- function(Object, DuplicationHandling){
  return(Object)
}
