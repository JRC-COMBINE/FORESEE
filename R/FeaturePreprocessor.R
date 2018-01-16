#' Preprocess the Inputs of both TrainObject and TestObject for Modeling Drug Response
#'
#' The FeaturePreprocessor converts the original features into predictive features that are defined by FeaturePreprocessing.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param TestObject Object that contains all data that the model is to be tested on, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param FeaturePreprocessing Method for preprocessing the inputs of the model, such as z-score, principal component analysis, PhysioSpace similarity, etc.

#' @return \item{TrainObject}{The TrainObject with preprocessed features.}
#'         \item{TestObject}{The TestObject with preprocessed features.}
#' @export

FeaturePreprocessor <- function(TrainObject, TestObject, FeaturePreprocessing){
  return(TrainObject, TestObject)
}
