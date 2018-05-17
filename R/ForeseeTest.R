#'  Test a Drug Efficacy Prediction Model on a TestObject
#'
#' ForeseeTest applies the machine learning based model ForeseeModel that has been trained on the data of a FORESEE TrainObject to a FORESEE TestObject to evaluate the Predictability of Drug Efficacy.
#' First, the Foreseer applies the ForeseeModel to the test data to gain the predicted response Foreseen.
#' Second, the Validator evaluates the performance of the model by comparing the predicted response Foreseen to the reported, true response.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type) and drug response data
#' @param ForeseeModel Model that has been trained on a TrainObject with the function ForeseeTrain.
#' @param Evaluation Measure for evaluating the model performance, such as ROC-Curve, AUC or p-value of ROC-Curve, Rsquared, MSE, Correlation, F-Test, etc.
#' Get all possible values with listInputOptions("Validator").
#' Instead of choosing one of the implemented options, a user-defined function can be used as an input.
#' @param BlackBox BlackBox used for training ForeseeModel.
#' @return \item{Performance}{Evaluation Measure of the Predictability of the ForeseeModel trained on the TrainObject and tested on the TestObject.}
#'         \item{Foreseen}{Predicted drug response of the TestObject obtained by applying the ForeseeModel to the molecular data of the TestObject.}
#' @export

#########################
# This file is part of the FORESEE R-package
# File authors: Lisa-Katrin Turnhoff <turnhoff@combine.rwth-aachen.de> and Ali Hadizadeh Esfahani <hadizadeh@combine.rwth-aachen.de>
# Distributed under the GNU General Public License v3.0.(http://www.gnu.org/licenses/gpl-3.0.html)
#########################

ForeseeTest <- function(TestObject, ForeseeModel, Evaluation = "rocauc", BlackBox = "ridge"){


    #################################################################################################################################
    # 1. Applying the model to test data

    Foreseen <- Foreseer(TestObject, ForeseeModel, BlackBox)
    #################################################################################################################################

    # 2. Calculating Performance

    Performance <- Validator(Foreseen, TestObject, Evaluation)

    #################################################################################################################################

    assign("Performance", value = Performance, envir = parent.frame())
    assign("Foreseen", value = Foreseen, envir = parent.frame())


}

