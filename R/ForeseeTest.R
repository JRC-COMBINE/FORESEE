#'  Test a Drug Efficacy Prediction Model on a TestObject
#'
#' ForeseeTest applies the machine learning based model ForeseeModel that has been trained on the data of a FORESEE TrainObject to a FORESEE TestObject to evaluate the Predictability of Drug Efficacy.
#' First, the Foreseer applies the ForeseeModel to the test data to gain the predicted response Foreseen.
#' Second, the Validator evaluates the performance of the model by comparing the predicted response Foreseen to the reported, true response.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param ForeseeModel Model that has been trained on a TrainObject with ForeseeTrain.
#' @param Evaluation Measure for evaluating the model performance, such as ROC-Curve, AUC or p-value of ROC-Curve, Rsquared, MSE, Correlation, F-Test, etc.
#' @param BlackBox Modeling algorithm for training, such as linear regression, elastic net, lasso regression, ridge regression, tandem, support vector machines, random forests, user defined functions, etc.

#' @return \item{Performance}{Evaluation Measure of the Predictability of the ForeseeModel trained on the TrainObject and tested on the TestObject.}
#'         \item{Foreseen}{Predicted drug response of the TestObject obtained by applying the ForeseeModel.}
#' @export

ForeseeTest <- function(TestObject, ForeseeModel, BlackBox, Evaluation){

  TestObject_test<- as.data.frame(as.matrix(t(TestObject$GeneExpression)))
  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

  Foreseen <- predict(object = ForeseeModel, newdata = TestObject_test)
  # Update Objects in the Environment
  Foreseen <<- Foreseen
  Performance <<- 0
  # Plots
}

