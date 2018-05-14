#' Genomics of Drug Sensitivity in Cancer or GDSC
#'
#' A dataset containing drug responses on cancer celllines.
#'
#'
#' @format ForeseeCell object
#' \describe{
#'   \item{GeneExpression}{Gene Expression Matrix}
#'   \item{AUC}{...}
#'   ...
#' }
#' @source \url{http://www.cancerrxgene.org}
"GDSC"


#' #' Landmark Genes of the Broad Institute (LINCS)
#' #'
#' #' A list of a subset of informative genes (ENTREZ IDs) determined by the Broad Institute.
#' #' Landmark genes were selected as those widely expressed across lineage and were found to have good predictive power for inferring the expression of other genes that are not directly measured in the L1000 assay.
#' #'
#' #' @format integers
#' #' @source \url{https://clue.io/command?q=/gene-space%20lm
#' #' }
#' "LM_genes_entrez"
#'
#'
#'
#' #' Human Housekeeping Genes
#' #'
#' #' A list of genes (ENTREZ IDs) that are involved in basic cell maintenance and, therefore, are expected to maintain constant expression levels in all cells and conditions.
#' #' Published by E. Eisenberg and E.Y. Levanon, Trends in Genetics, 29 (2013) "Human housekeeping genes revisited"
#' #'
#' #' @format integers
#' #' @source \url{https://www.tau.ac.il/~elieis/HKG/}
#' "HK_genes_entrez"


