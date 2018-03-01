#' Transform Drug Response Data
#'
#' The CellResponseProcessor transforms the response data of the TrainObject for prediction.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param DrugName Name of the drug whose efficacy is supposed to be predicted with the model
#' @param CellResponseType Format of the drug response data of the TrainObject, such as IC50, AUC, GI50, etc., that is used for prediction
#' @param CellResponseTransformation Method that is to be used to transform the drug response data of the TrainObject:
#' the function 'powertransform' power transforms the drug response data,
#' the function 'logarithm' returns the natural logarithm of the drug response data,
#' the function 'binarization_kmeans' returns a binarized drug response vector based on 2 kmeans clusters,
#' the function 'binarization_cutoff' returns a binarized drug response vector based on a cutoff at the median,
#' the function 'none' returns the unchanged drug response data,
#' the function 'user-defined function' is determined by the function in the input
#' @return \item{TrainObject}{The TrainObject with preprocessed drug response data.}
#' @export


# ToDo
# User-defined input not implemented yet (see commented, not-working part CellResponseProcessor.function)
# Include cutoff binarization for transformation type (remove package, simply divide input vector at median)
# Inlude options for other response types (what if user choses values other than AUC or IC50?)

CellResponseProcessor <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){
  UseMethod("CellResponseProcessor", object = CellResponseTransformation)
  return(TrainObject)
}

CellResponseProcessor.character <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){
  class(CellResponseTransformation) <- CellResponseTransformation;
  UseMethod("CellResponseProcessor", object = CellResponseTransformation)
}

### NOT WORKING YET
# CellResponseProcessor.function <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){
#   UseMethod("CellResponseProcessor", object = CellResponseTransformation)
#   message("The used-defined function is applied")
#
#   Object_withDrugResponse <- GetCellResponseData(TrainObject = TrainObject, DrugName = DrugName, CellResponseType = CellResponseType)
#   # ...
#   # Update TrainObject in the Environment
#   TrainObject <<- Object_withDrugResponse
# }



################################################################################
### Function "powertransform" to powertransform the chosen drug response data
CellResponseProcessor.powertransform <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){

  # Load Package for Power Transform
  require(car)

  # Number of Cell Lines before adjusting with drug response data
  dim_before <- dim(TrainObject$GeneExpression)[2]

  # Extract drug response of interest
  Object_withDrugResponse <- GetCellResponseData(TrainObject = TrainObject, DrugName = DrugName, CellResponseType = CellResponseType)

  # Do powertransform of drug response data

  # Powertransform needs all inputs to be positive
    if(min(Object_withDrugResponse$DrugResponse, na.rm = TRUE) < 0) {
      offset <- -min(Object_withDrugResponse$DrugResponse, na.rm = TRUE) + 1
      Object_withDrugResponse$DrugResponse <- Object_withDrugResponse$DrugResponse + offset
    }

  TransForm <- powerTransform(Object_withDrugResponse$DrugResponse)$lambda
  Object_withDrugResponse$DrugResponse <- Object_withDrugResponse$DrugResponse^TransForm

  # Number of Cell Lines after adjusting with drug response data
  dim_after <- dim(Object_withDrugResponse$GeneExpression)[2]

  # Prints the reduction of gene names
  message(paste0("CellResposeProcessor added the new matrix 'Drug Response' to the ForeseeCell Object, which includes power transformed ",CellResponseType," response information about ",DrugName,"."))
  message(paste0("The number of cell lines in the ForeseeCell Object was reduced from ",dim_before," to ",dim_after,"."))

  # Returns the new TrainObj
  return(Object_withDrugResponse)
}


################################################################################
### Function "logarithm" to logarithm the chosen drug response data
CellResponseProcessor.logarithm <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){

  # Number of Cell Lines before adjusting with drug response data
  dim_before <- dim(TrainObject$GeneExpression)[2]

  # Extract drug response of interest
  Object_withDrugResponse <- GetCellResponseData(TrainObject = TrainObject, DrugName = DrugName, CellResponseType = CellResponseType)


  # Log needs all inputs to be positive, otherwise only NAs are returned
  if(min(Object_withDrugResponse$DrugResponse, na.rm = TRUE) < 0) {
    offset <- -min(Object_withDrugResponse$DrugResponse, na.rm = TRUE) + 1
    Object_withDrugResponse$DrugResponse <- Object_withDrugResponse$DrugResponse + offset
  }

  # Do logarithm of drug response data
  Object_withDrugResponse$DrugResponse <- log(Object_withDrugResponse$DrugResponse)

  # Number of Cell Lines after adjusting with drug response data
  dim_after <- dim(Object_withDrugResponse$GeneExpression)[2]

  # Prints the reduction of gene names
  message(paste0("CellResposeProcessor added the new matrix 'Drug Response' to the ForeseeCell Object, which includes natural logarithmic ",CellResponseType," response information about ",DrugName,"."))
  message(paste0("The number of cell lines in the ForeseeCell Object was reduced from ",dim_before," to ",dim_after,"."))

  # Returns the new TrainObj
  return(Object_withDrugResponse)
}


################################################################################
### Function "binarization_kmeans" to binarization the chosen drug response data
### Uses the kmeans algorithm of the package Binarize to find two clusters in the data

CellResponseProcessor.binarization_kmeans <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){

  require(Binarize)

  # Number of Cell Lines before adjusting with drug response data
  dim_before <- dim(TrainObject$GeneExpression)[2]

  # Extract drug response of interest
  Object_withDrugResponse <- GetCellResponseData(TrainObject = TrainObject, DrugName = DrugName, CellResponseType = CellResponseType)

  # Do kmeans binarization of drug response data
  Object_withDrugResponse$DrugResponse <- binarize.kMeans(Object_withDrugResponse$DrugResponse)@binarizedMeasurements


  # Number of Cell Lines after adjusting with drug response data
  dim_after <- dim(Object_withDrugResponse$GeneExpression)[2]

  # Prints the reduction of gene names
  message(paste0("CellResposeProcessor added the new matrix 'Drug Response' to the ForeseeCell Object, which includes binarized ",CellResponseType," response information about ",DrugName,"."))
  message(paste0("The number of cell lines in the ForeseeCell Object was reduced from ",dim_before," to ",dim_after,"."))

  # Returns the new TrainObj
  return(Object_withDrugResponse)
}

################################################################################
### Function "binarization_cutoff" to binarization the chosen drug response data
### Uses the kmeans algorithm of the package Binarize to find two clusters in the data

CellResponseProcessor.binarization_cutoff <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){

  require(bootnet)

  # Number of Cell Lines before adjusting with drug response data
  dim_before <- dim(TrainObject$GeneExpression)[2]

  # Extract drug response of interest
  Object_withDrugResponse <- GetCellResponseData(TrainObject = TrainObject, DrugName = DrugName, CellResponseType = CellResponseType)

  # Do binarization of drug response data with median as cutoff
  Object_withDrugResponse$DrugResponse <- binarize(x=Object_withDrugResponse$DrugResponse, split = "median", removeNArows = TRUE)$x


  # Number of Cell Lines after adjusting with drug response data
  dim_after <- dim(Object_withDrugResponse$GeneExpression)[2]

  # Prints the reduction of gene names
  message(paste0("CellResposeProcessor added the new matrix 'Drug Response' to the ForeseeCell Object, which includes binarized ",CellResponseType," response information about ",DrugName,"."))
  message(paste0("The number of cell lines in the ForeseeCell Object was reduced from ",dim_before," to ",dim_after,"."))

  # Returns the new TrainObj
  return(Object_withDrugResponse)
}


################################################################################
### Function "none" to use the raw drug response data
CellResponseProcessor.none <- function(TrainObject, DrugName, CellResponseType, CellResponseTransformation){

  # Number of Cell Lines before adjusting with drug response data
  dim_before <- dim(TrainObject$GeneExpression)[2]

  # Extract drug response of interest
  Object_withDrugResponse <- GetCellResponseData(TrainObject = TrainObject, DrugName = DrugName, CellResponseType = CellResponseType)

  # Don't do anything to drug response data

  # Number of Cell Lines after adjusting with drug response data
  dim_after <- dim(Object_withDrugResponse$GeneExpression)[2]

  # Prints the reduction of gene names
  message(paste0("CellResposeProcessor added the new matrix 'Drug Response' to the ForeseeCell Object, which includes binarized ",CellResponseType," response information about ",DrugName,"."))
  message(paste0("The number of cell lines in the ForeseeCell Object was reduced from ",dim_before," to ",dim_after,"."))

  # Returns the new TrainObj
  return(Object_withDrugResponse)
}
