## Script that will download, process and add GSE6434 to this package
# by esfahani@aices.rwth-aachen.de, 22.1.2018

require(affy)
GSE6434 <- justRMA(celfile.path = "./data-raw/GSE6434_RAW/")

GSE6434GEX <- exprs(GSE6434)

require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("affy_hg_u95av2","entrezgene"),
                     filters = "affy_hg_u95av2", values = rownames(GSE6434GEX), mart = humaRt)

rownames(GSE6434GEX) <- ConvTabelle$entrezgene[match(rownames(GSE6434GEX),ConvTabelle$affy_hg_u95av2)]
GSE6434GEX <- GSE6434GEX[!is.na(rownames(GSE6434GEX)),] #1636 probes lost in Entrez conversion

ActualAnnotation <- character(length = 24)
ActualAnnotation[c(1,2,4,5,6,9,10,11,12,13,16,18,22,24)] <- "Resistant"
ActualAnnotation[c(3,7,8,14,15,17,19,20,21,23)] <- "Sensitive"

Annotation <- logical(length = 24)
Annotation[c(1,2,4,5,6,9,10,11,12,13,16,18,22,24)] <- FALSE
Annotation[c(3,7,8,14,15,17,19,20,21,23)] <- TRUE

names(Annotation) <- ActualAnnotation


#Making of the ForeseePatient object:
GSE6434expressionSet <- GSE6434 #So I wouldn't overwrite it
GSE6434 <- list()
class(GSE6434) <- "ForeseePatient"
GSE6434[["GeneExpression"]] <- GSE6434GEX
GSE6434[["Annotation"]] <- Annotation

devtools::use_data(GSE6434, overwrite = F)






