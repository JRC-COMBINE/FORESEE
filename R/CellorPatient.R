#' Tests if the Object is a Cell or a Patient Based on the Drug Response Type
#'
#' The CellorPatient Function tests if the object is a cell or a patient based on the type of the drug response annotation.
#'
#' @param Object FORESEE TrainObject or TestObject that should be analyzed.
#'
#' @return \item{CellorPatient_status}{The Status of the FORESEE Object that states, whether it is a cell or a patient.}
#' @export

CellorPatient <- function(Object){
  if(any(class(Object)==c("ForeseeCell","ForeseePatient"))) {
    warning(paste("Input of CellorPatient function should be a Foresee Object! returning NULL"))
    return(NULL)
  }
  return(ifelse(test = class(GSE6434)=="ForeseeCell", yes = T, no = F))
}
