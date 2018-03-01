#' List Input Options for FORESEE Methods
#'
#' listInputOptions returns possible input values in ForeseeTrain (CellResponseTransformation, DuplicationHandling,
#' HomogenizationMethod, GeneFilter, FeaturePreprocessing and BlackBox) and in ForeseeTest (Evaluation).
#'
#' @param METHOD Character string of the function corrisponding to the input of ForeseeTrain or ForeseeTest.
#' Check the help of Foreseetrain and ForeseeTest for more info.
#'
#' @return \item{Character vector of all possible inputs}
#' @examples
#' listInputOptions("CellResponseProcessor")
#' listInputOptions("DuplicationHandler")
#' listInputOptions("Homogenizer")
#' listInputOptions("Validator")
#' @export

listInputOptions <- function(METHOD){
  listOfPossibilities <- methods(METHOD)
  listOfPossibilitiesWithOutMETHODName <- sapply(strsplit(listOfPossibilities, split = ".", fixed = T), function(x) x[2])
  listOfPossibilitiesWithOutMETHODName <- setdiff(listOfPossibilitiesWithOutMETHODName, c("default", "character")) # "default" and "character" are used in implementation but nor actually input options
  if(any(listOfPossibilitiesWithOutMETHODName == "function")){
    listOfPossibilitiesWithOutMETHODName[listOfPossibilitiesWithOutMETHODName=="function"] <- "User Defined Function"
  }
  return(listOfPossibilitiesWithOutMETHODName)
}
