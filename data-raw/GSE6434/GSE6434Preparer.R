## Script that will download, process and add GSE6434 to FORESEE package
# by esfahani@aices.rwth-aachen.de, 22.1.2018


##Raw CEL files of GSE6434 must be downloaded from Gene Expression Omnibus (GEO) via the
#link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6434, and untared in GSE6434-raw folder.
##You can do this process automatically by running the next lines:
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE6434&format=file",
              destfile = "./data-raw/GSE6434/GSE6434_RAW.tar") #Downloading the CEL files.
untar("./data-raw/GSE6434/GSE6434_RAW.tar", exdir = "./data-raw/GSE6434/GSE6434-raw/") #Extract files
unlink("./data-raw/GSE6434/GSE6434_RAW.tar") #Removing the tar file


###Preparing Gene Expression Data:
##We use affy package to do RMA of the CEL files:
require(affy)
GSE6434 <- justRMA(celfile.path = "./data-raw/GSE6434/GSE6434-raw/")

#Extracting gene expression matrix:
GSE6434GEX <- exprs(GSE6434)

#Converting affy gene IDs to Entrez IDs:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("affy_hg_u95av2","entrezgene"),
                     filters = "affy_hg_u95av2", values = rownames(GSE6434GEX), mart = humaRt)

rownames(GSE6434GEX) <- ConvTabelle$entrezgene[match(rownames(GSE6434GEX),ConvTabelle$affy_hg_u95av2)]
GSE6434GEX <- GSE6434GEX[!is.na(rownames(GSE6434GEX)),]


###Making the Annotation:
#We extracted the annotation based on information on https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6434:
ActualAnnotation <- character(length = 24)
ActualAnnotation[c(1,2,4,5,6,9,10,11,12,13,16,18,22,24)] <- "Resistant"
ActualAnnotation[c(3,7,8,14,15,17,19,20,21,23)] <- "Sensitive"

#Since basically we have a binary value for annotaion, we make annotation a logical vector, with the actual
#meanings of True and False in 'names' of the logical vector:
Annotation <- logical(length = 24)
Annotation[c(1,2,4,5,6,9,10,11,12,13,16,18,22,24)] <- FALSE
Annotation[c(3,7,8,14,15,17,19,20,21,23)] <- TRUE
names(Annotation) <- ActualAnnotation


###Making and saving the Foresee object:
GSE6434 <- list()
class(GSE6434) <- "ForeseePatient"
GSE6434[["GeneExpression"]] <- GSE6434GEX
GSE6434[["Annotation"]] <- Annotation

devtools::use_data(GSE6434, overwrite = F)






