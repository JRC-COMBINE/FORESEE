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


    #################################################################################################################################
    # 1. Applying the model to test data

    Foreseen <- Foreseer(TestObject, ForeseeModel)
    #################################################################################################################################


    # 1. Applying the model to test data

    # Validation <- Validator(Foreseen,...  )
    #################################################################################################################################




}

