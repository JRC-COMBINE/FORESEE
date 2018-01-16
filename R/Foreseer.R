#' Apply the trained ForeseeModel on a TestObject
#'
#' The Foreseer applies the ForeseeModel that was trained on a FORESEE TrainObject to the test data to gain a prediction of the TestObject's response.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param ForeseeModel Model that has been trained on a TrainObject with ForeseeTrain.
#'
#' @return \item{Foreseen}{Predicted drug response of the TestObject obtained by applying the ForeseeModel.}
#' @export

Foreseer <- function(TestObject, ForeseeModel){
  Foreseen<-0
  return(Foreseen)
}


