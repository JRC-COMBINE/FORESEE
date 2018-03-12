#' List Input Options for FORESEE Methods
#'
#' listInputOptions returns possible input values in ForeseeTrain (CellResponseTransformation, DuplicationHandling,
#' HomogenizationMethod, GeneFilter, FeaturePreprocessing and BlackBox) and in ForeseeTest (Evaluation).
#'
#' @param METHOD Character string of the function corrisponding to the input of ForeseeTrain or ForeseeTest.
#' Check the help of Foreseetrain and ForeseeTest for more info.
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

#' List All Cell lines Inside a ForeseeTrain Object
#'
#' listCellLines returns all cell lines available in a ForeseeTrain Object.
#'
#' @param OBJ A ForeseeTrain object that you want to extract its containing cell lines.
#' @return \item{Character vector of all cell line names in OBJ.}
#' @examples
#' listCellLines(GDSC)
#' @export

listCellLines <- function(OBJ){
  return(colnames(OBJ$GeneExpression))
}

#' List All Drugs Inside a ForeseeTrain Object
#'
#' listDrugs returns all possible drugs that were tested and has reponses available in a ForeseeTrain Object.
#'
#' @param OBJ A ForeseeTrain object that you want to extract all possible drugs it has response information on.
#' @return \item{Character vector of all drug names in OBJ.}
#' @examples
#' listDrugs(GDSC)
#' @export

listDrugs <- function(OBJ){
  if(!any(names(OBJ) == "IC50")) stop("listDrugs is not implemented completely!!")
  return(colnames(OBJ$IC50)) ## This will break when OBJ doesn't have IC50 -> ToDO: have to have a better 'Universal' slot available in all Train Objects
}
