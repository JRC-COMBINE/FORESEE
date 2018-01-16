#' Transform Drug Response Data
#'
#' The CellResponseProcessor transforms the response data of the TrainObject for prediction.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param DrugName Name of the drug whose efficacy is supposed to be predicted with the model
#' @param CellResponseType Format of the drug response data of the TrainObject, such as IC50, AUC, GI50, etc., that is used for prediction
#' @param CellResponseTransformation Method that is to be used to transform the drug response data of the TrainObject, such as power transform, logarithm, binarization, user defined functions, etc.
#'
#' @return \item{TrainObject}{The TrainObject with preprocessed drug response data.}
#' @export

CellResponseProcessor <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){
  return(TrainObject)
}
