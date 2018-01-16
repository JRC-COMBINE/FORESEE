#' Homogenize a Pair of two FORESEE Objects
#'
#' The Homogenizer function removes batch effects and homogenizes the data of a given TrainObject and TestObject.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param TestObject Object that contains all data that the model is to be tested on, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param HomogenizationMethod Method for homogenizing data of the TrainObject and TestObject, such as ComBat, quantile normalization, limma, RUV, etc.
#'
#' @return \item{TrainObject}{The TrainObject with homogenized features according to the chosen TestObject.}
#'         \item{TestObject}{The TestObject with homogenized features according to the chosen TrainObject.}
#' @export

Homogenizer <- function(TrainObject, TestObject, HomogenizationMethod){
  return(TrainObject, TestObject)
}
