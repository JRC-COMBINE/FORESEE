#' Get Cell Line Drug Response from FORESEE Object
#'
#' The GetCellResponseData Function extracts the features of the drug response data that are relevant to predict the drug response.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param DrugName Name of the drug whose efficacy is supposed to be predicted with the model
#' @param CellResponseType Format of the drug response data of the TrainObject, such as IC50, AUC, GI50, etc., that is used for prediction
#'
#' @return \item{TrainObject}{The TrainObject with extracted drug response data.}
#' @export

GetCellResponseData <- function(TrainObject, DrugName, CellResponseType){

  DrugResponse <- TrainObject[[CellResponseType]][,grep(DrugName,colnames(TrainObject[[CellResponseType]]))]
  DrugResponse <- DrugResponse[is.na(DrugResponse)==FALSE]
  TrainObject[["DrugResponse"]] <- DrugResponse

  return(TrainObject)
}
