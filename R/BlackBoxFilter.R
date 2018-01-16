#' Train a Black Box Model for Drug Efficacy Prediction
#'
#' The BlackBoxFilter applies a machine learning algorithm to the data of the TrainObject to create a model that is predictive of the drug response.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param BlackBox Modeling algorithm for training, such as linear regression, elastic net, lasso regression, ridge regression, tandem, support vector machines, random forests, user defined functions, etc.
#' @param nfoldCrossvalidation # folds to use for crossvalidation while training the model. If put to zero, the complete data of the TrainObject is used for training.

#' @return \item{ForeseeModel}{A black box model trained on the TrainObject data that can be applied to new test data.}
#'         \item{TrainObject}{The TrainObject that was used to train the model.}
#' @export

BlackBoxFilter <- function(TrainObject, BlackBox, nfoldCrossvalidation){
  ForeseeModel<-0
  return(ForeseeModel,TrainObject)
}
