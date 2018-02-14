#' Validate the Performance of a ForeseeModel on a new TestObject
#'
#' The Validator evaluates the performance of the model by comparing the predicted response Foreseen to the reported, true response from the TestObject's Annotation.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, especially the true, measured drug response.
#' @param Foreseen Predicted drug response of the TestObject obtained by applying the ForeseeModel.
#' @param Evaluation Measure for evaluating the model performance, such as ROC-Curve, AUC or p-value of ROC-Curve, Rsquared, MSE, Correlation, F-Test, etc.

#' @return \item{Performance}{Evaluation Measure of the Predictability of the ForeseeModel trained on the TrainObject and tested on the TestObject.}
#' @export

Validator <- function(Foreseen, TestObject, Evaluation){
  UseMethod("Validator", object = HomogenizationMethod)
}

Validator.character <- function(Foreseen, TestObject, Evaluation){
  class(Evaluation) <- Evaluation;
  UseMethod("Validator", object = Evaluation)
  return(Performance)
}

################################################################################
### Function "function" applies the function in "Evaluation" to predicted values
Validator.function <- function(Foreseen, TestObject, Evaluation) {

  return(Performance)
}

################################################################################
### Function "function" applies the function in "Evaluation" to predicted values
Validator.rocauc <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <-TestObject$Annotation
  if(is.numeric(ANNOTATIONS)){
    message("Annotation of the test set is binarized for calculating ROC")
    ANNOTATIONS <- ifelse(ANNOTATIONS-median(ANNOTATIONS) > 0, TRUE, FALSE)
  }
  require(pROC)
  ROCObj <- roc(-(as.numeric(ANNOTATIONS)), Foreseen, direction="<")
  AUCofROC <- pROC::auc(pROC::roc(ANNOTATIONS, Foreseen))[[1]]
  return(AUCofROC)

  ### Plot ROC Curve
  jpeg(filename = filename_roc, width=10, height=10, units="in", res=600)
  plot(ROCObj, main = "Prediction of Patient Response after Training on Cell Line Data",legacy.axes = TRUE, print.auc = TRUE)
  dev.off()


}





