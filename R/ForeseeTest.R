#'  Test a Drug Efficacy Prediction Model on a TestObject
#'
#' ForeseeTest applies the machine learning based model ForeseeModel that has been trained on the features of a FORESEE TrainObject to a FORESEE TestObject to evaluate the Predictability of Drug Efficacy.
#' First, the Foreseer applies the ForeseeModel to the test data to gain the predicted response 'Foreseen'.
#' Second, the Validator evaluates the performance of the model by comparing the predicted response 'Foreseen' to the reported, true response.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc. ) and drug response data
#' @param ForeseeModel Model that has been trained on a TrainObject with the function ForeseeTrain.
#' @param Evaluation Measure for evaluating the model performance.
#' The option 'fpvalue' calculates the p value of an F test on a linear model between predictions and the actual annotations,
#' The option 'mse' calculates the mean square error of a linear model between predictions and the actual annotations,
#' The option 'pearson' calculates the pearson correlation between predictions and the actual annotations,
#' The option 'prauc' calculates the AUC of the precision-recall curve
#' The option 'rocauc' calculates the AUC of the ROC curve
#' The option 'rocpvalue' calculates the t.test p value of the ROC curve AUC versus 10000 permutated annotation values,
#' The option 'rsquared' calculates the fraction of variance explained by a linear model between predictions and actual annotations,
#' The option 'rsquared_adjusted' calculates the fraction of variance explained by a linear model between predictions and actual annotations, corrected the p-value of F-test,
#' The option 'spearman'. calculates the spearman correlation between predictions and the actual annotations.
#' The function 'listInputOptions("Validator")' returns a list of the possible options.
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

