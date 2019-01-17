## Script that will download, process and add GSE51373 to FORESEE package
# by esfahani@aices.rwth-aachen.de, 14.11.2018


##Raw CEL files of GSE51373 must be downloaded from Gene Expression Omnibus (GEO) via the
#link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51373, and untared in GSE51373-raw folder.
##You can do this process automatically by running the next lines:
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE51373&format=file",
              destfile = "data-raw/GSE51373/GSE51373-raw/GSE51373_RAW.tar") #Downloading the CEL files.
untar("data-raw/GSE51373/GSE51373-raw/GSE51373_RAW.tar", exdir = "data-raw/GSE51373/GSE51373-raw/") #Extract files
unlink("data-raw/GSE51373/GSE51373-raw/GSE51373_RAW.tar") #Removing the tar file


###Preparing Gene Expression Data:
#We use affy package to do RMA of the CEL files:
require(affy)
GSE51373 <- justRMA(celfile.path = "data-raw/GSE51373/GSE51373-raw/")
#We also download GSE51373 object (already normalized) from GEO to use its anotations:
library(GEOquery)
GSE51373PreRMAed <- getGEO(GEO = "GSE51373", getGPL = FALSE)#It's MAS5!
GSE51373PreRMAed <- GSE51373PreRMAed$GSE51373_series_matrix.txt.gz#just want to use it's annots
#Extracting gene expression matrix:
GSE51373GEX <- exprs(GSE51373)

#Converting affy gene IDs to Entrez IDs (using GPL object from GEO):
GPL570 <- getGEO(GEO = "GPL570")
rownames(GSE51373GEX) <- GPL570@dataTable@table$ENTREZ_GENE_ID[match(rownames(GSE51373GEX),GPL570@dataTable@table$ID)]
rownames(GSE51373GEX) <-
  sapply(strsplit(rownames(GSE51373GEX), split = " /// ", fixed = TRUE), function(x)
    x[1]) #In multiple matches, we get the first Entrez ID.


###Making the Annotation:
#We extracted the annotation based on information on the original publication of
#this data set (Koti et. al., BMC Cancer 2013,
#https://bmccancer.biomedcentral.com/articles/10.1186/1471-2407-13-549):
ActualAnnotation <- character(length = 28)
ActualAnnotation[c(4,6,3,20,2,16,17,15,18,28)] <- "Resistant"
ActualAnnotation[c(21,9)] <- "Partially Resistant"
ActualAnnotation[c(7,5,8,19,11,12,1,22,14,10,23,13,24,25,27,26)] <- "Sensitive"

#Since basically we have a binary value for annotaion, we make annotation a logical vector, with the actual
#meanings of True and False in 'names' of the logical vector:
Annotation <- logical(length = 28)
Annotation[c(4,6,3,20,2,16,17,15,18,28,21,9)] <- FALSE
Annotation[c(7,5,8,19,11,12,1,22,14,10,23,13,24,25,27,26)] <- TRUE
names(Annotation) <- ActualAnnotation


###Making and saving the Foresee object:
GSE51373 <- list()
class(GSE51373) <- "ForeseePatient"
GSE51373[["GeneExpression"]] <- GSE51373GEX
GSE51373[["Annotation"]] <- Annotation
GSE51373[["ExtraAnnotation"]] <- GSE51373PreRMAed@phenoData@data

devtools::use_data(GSE51373, overwrite = F)

