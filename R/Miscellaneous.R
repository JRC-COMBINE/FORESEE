#' List Input Options for FORESEE Methods
#'
#' listInputOptions returns possible input arguments in ForeseeTrain (for DrugName, CellResponseType, InputDataTypes, TrainingTissue,
#' TestingTissue, CellResponseTransformation, DuplicationHandling, HomogenizationMethod, GeneFilter, FeaturePreprocessing
#' and BlackBox arguments) and in ForeseeTest (for Evaluation argument).
#'
#' @param FunArgument Character string of the input argument (for DrugName, CellResponseType, InputDataTypes, TrainingTissue
#' and TestingTissue) or the function corresponding to the input of ForeseeTrain (CellResponseTransformation, DuplicationHandling,
#' HomogenizationMethod, GeneFilter, FeaturePreprocessing and BlackBox) or ForeseeTest (Evaluation).

#' Check the help of ForeseeTrain and ForeseeTest for more information.
#' @param OBJ A ForeseeTrain or ForeseeTest object that you want to extract options for. Only necessary for when FunArgument is
#' DrugName, CellResponseType, InputDataTypes, TrainingTissue or TestingTissue.
#' @return \item{Character vector of all possible inputs}{}
#' @examples
#' listInputOptions("CellResponseProcessor")
#' listInputOptions("DuplicationHandler")
#' listInputOptions("Homogenizer")
#' listInputOptions("Validator")
#' listInputOptions("DrugName", CCLE)[1:10]
#' listInputOptions("InputDataTypes", GAO)
#' listInputOptions("CellResponseType", GDSC)
#' listInputOptions("TrainingTissue", CCLE)
#' listInputOptions("TestingTissue", GSE33072_erlotinib)
#' @export

listInputOptions <- function(FunArgument, OBJ){
  if(FunArgument=="DrugName"){
    return(colnames(OBJ[[as.character(OBJ$ResponseTypes$Name[1])]]))
  } else if(FunArgument=="CellResponseType"){
    return(as.character(OBJ$ResponseTypes$Name))
  } else if(FunArgument=="InputDataTypes"){
    return(as.character(OBJ$InputTypes$Name))
  } else if(FunArgument=="TrainingTissue"){
    if(any(names(OBJ) == "TissueInfo") & any(names(OBJ$TissueInfo) == "Site")){
      return(as.character(unique(OBJ$TissueInfo$Site)))
    } else {
      warning("Object doesn't have tissue site information, the only acceptable option for TrainingTissue is 'all' for using all samples.")
      return("all")
    }
  } else if(FunArgument=="TestingTissue"){
    if(any(names(OBJ) == "TissueInfo") & any(names(OBJ$TissueInfo) == "Site")){
      return(as.character(unique(OBJ$TissueInfo$Site)))
    } else {
      warning("Object doesn't have tissue site information, the only acceptable option for TestingTissue is 'all' for using all samples.")
      return("all")
    }
  } else {
    listOfPossibilities <- methods(FunArgument)
    listOfPossibilitiesWithOutMETHODName <- sapply(strsplit(listOfPossibilities, split = ".", fixed = T), function(x) x[2])
    listOfPossibilitiesWithOutMETHODName <- setdiff(listOfPossibilitiesWithOutMETHODName, c("default", "character")) # "default" and "character" are used in implementation but nor actually input options
    if(any(listOfPossibilitiesWithOutMETHODName == "function")){
      listOfPossibilitiesWithOutMETHODName[listOfPossibilitiesWithOutMETHODName=="function"] <- "User Defined Function"
    }
    return(listOfPossibilitiesWithOutMETHODName)
  }
}

#' List All Cell lines Inside a ForeseeTrain Object
#'
#' listCellLines returns all cell lines (or sample names in case of a xenograft data set)
#' available in a ForeseeTrain Object (All cell lines that are included in gene expression matrix
#' of ForeseeTrain Object).
#'
#'
#' @param OBJ A ForeseeTrain object of which you want to extract its containing cell lines.
#' @return \item{Character vector of all cell line names in OBJ.}{}
#' @examples
#' listCellLines(GDSC)
#' listCellLines(WITKIEWICZ)
#' @export

listCellLines <- function(OBJ){
  return(colnames(OBJ$GeneExpression))
}

#' List All Drugs Inside a ForeseeTrain Object
#'
#' listDrugs returns all possible drugs that were tested and have reponses available in a ForeseeTrain Object.
#'
#' @param OBJ A ForeseeTrain object that you want to extract all possible drugs it has response information on.
#' @return \item{Character vector of all drug names in OBJ.}{}
#' @examples
#' listDrugs(GDSC)
#' listDrugs(GAO)
#' @export

listDrugs <- function(OBJ){
  if(class(OBJ) != "ForeseeCell") stop(paste("'listDrugs' is only applicable for ForeseeCell objects!, class of input OBJ is",class(OBJ),"!"))
  return(listInputOptions("DrugName",OBJ))
}

#' Checking CellResponseType Availability
#'
#' Checking If CellResponseType is Available in Object
#'
#' @param OBJ A ForeseeTrain object that you want to check the CellResponseType availability in.
#' @param RESP CellResponseType to be checked.
#' @return \item{Returns an invisible TRUE if RESP is available in OBJ, if not available an error will be generated}{}
#' @examples
#' CellResponseTypeAvailabilityCheck(DAEMEN,"GI50")
#' @export

CellResponseTypeAvailabilityCheck <- function(OBJ, RESP){
  if(!any(OBJ$ResponseTypes$Name == RESP)){
    stop(paste(RESP,"is not available as a CellResponseType in the TrainObject, You can check the names and descriptions of all available CellResponseTypes of your TrainObject in its 'ResponseTypes' component (e.g. by TrainObject$ResponseTypes)."))
  }
  invisible(TRUE)
}


#' Test if the Object is a Cell or a Patient Based on the Drug Response Type
#'
#' The CellorPatient Function tests if the object is a cell or a patient based on the type of the drug response annotation.
#'
#' @param Object FORESEE TrainObject or TestObject that should be analyzed.
#'
#' @return \item{CellorPatient_status}{The Status of the FORESEE Object that states, whether it is a cell (true) or a patient (false).}
#' @examples
#' CellorPatient(GDSC)
#' CellorPatient(GSE33072_sorafenib)
#' @export

CellorPatient <- function(Object){
  if(!any(class(Object)==c("ForeseeCell","ForeseePatient"))) {
    warning(paste("Input of CellorPatient function should be a Foresee Object! returning NULL"))
    return(NULL)
  }
  return(ifelse(test = class(Object)=="ForeseeCell", yes = T, no = F))
}


#' Loading/Attaching (and Installing) a Package
#'
#' requireForesee is the same as 'require' from the base package, except in case of a missing package, it
#' tries to install it via 'biocLite' from Bioconductor. For installation to work, R needs to have access
#' to the internet (more precisely "https://bioconductor.org/biocLite.R" should be accessable to R).
#'
#' @param package the name of a package to be loaded and attached (and installed).
#'
#' @return requireForesee returns (invisibly) a logical indicating whether the required package was available (before installation attempts).
#' @examples
#' requireForesee(ranger)
#' @export

requireForesee <- function(package){
  package <- as.character(substitute(package)) # In case user didn't provide a character but just typed the name of the package!
  suppressWarnings( # Cause 'require' makes a warning if package is not available
    loadSucceed <- require(package, quietly = TRUE, character.only = TRUE)
  )
  if(!loadSucceed){
    if(!interactive()){ #'interactive' cause there would be a user to ok the new package installation
      stop("R is not used interactively, 'requireForesee' needs to install a new package, which is only possible in interactive R!")
    } else {
      usrAns <- readline(prompt = paste("Packege",package,"is needed and not installed, would you like FORESEE to install it for you?(type y or yes): "))
      if(identical(usrAns, "y") | identical(usrAns, "yes")) {
        if(identical(package,"PhysioSpaceMethods")){
          requireForesee(devtools)
          install_github(repo = "JRC-COMBINE/PhysioSpaceMethods", build_vignettes = TRUE)
        } else {
          source("https://bioconductor.org/biocLite.R")
          biocLite(package)
        }
        library(package, character.only = TRUE)
      } else {
        stop(paste("Installation of",package,"was aborted by the user"))
      }
    }
  }
  invisible(loadSucceed)
}
