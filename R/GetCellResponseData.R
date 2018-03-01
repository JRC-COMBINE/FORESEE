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

  #Checking to see user's DrugName is available:
  if(!any(colnames(TrainObject[[CellResponseType]])==DrugName)){
    message(paste(DrugName,"wasn't found, trying to match case-insensitive..."))
    if(any(tolower(colnames(TrainObject[[CellResponseType]]))==tolower(DrugName))){
      if(sum(colnames(tolower(TrainObject[[CellResponseType]]))==tolower(DrugName)) > 1){
        stop("More than one drug were matched! Can not continue with more than one drug to access!")
      } else {
        DrugName <- colnames(TrainObject[[CellResponseType]])[colnames(tolower(TrainObject[[CellResponseType]]))==tolower(DrugName)]
      }
      message(paste(DrugName,"was found!"))
    } else {
      stop(paste("No drugs were found in TrainObject matching",DrugName,", Try listDrugs(TrainObject) to get all acceptable DrugNames in your TrainObject"))
    }
  }
  DrugResponse <- TrainObject[[CellResponseType]][,DrugName]
  DrugResponse <- DrugResponse[is.na(DrugResponse)==FALSE]
  TrainObject[["DrugResponse"]] <- DrugResponse

  # Return gene expression of only those cell lines that have drug response values
  CommonCelllines <- colnames(TrainObject$GeneExpression)[colnames(TrainObject$GeneExpression) %in% names(TrainObject$DrugResponse)]
  TrainObject$GeneExpression <- TrainObject$GeneExpression[,CommonCelllines]
  TrainObject$DrugResponse <- TrainObject$DrugResponse[CommonCelllines]

  return(TrainObject)
}
