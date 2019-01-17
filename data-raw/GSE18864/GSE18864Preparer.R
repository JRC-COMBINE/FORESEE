## Script that will download, process and add GSE18864 to FORESEE package
# by esfahani@aices.rwth-aachen.de, 27.2.2018


##Raw CEL files of GSE18864 must be downloaded from Gene Expression Omnibus (GEO) via the
#link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18864, (or from ArrayExpress
#using the this link: https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-18864/)
# and untared in GSE18864-raw folder.
##You can do this process automatically by running the next lines:
download.file(url = "https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-18864/E-GEOD-18864.raw.1.zip",
              destfile = "./data-raw/GSE18864/E-GEOD-18864.raw.1.zip") #Downloading the CEL files.
download.file(url = "https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-18864/E-GEOD-18864.raw.2.zip",
              destfile = "./data-raw/GSE18864/E-GEOD-18864.raw.2.zip") #Downloading the CEL files.
download.file(url = "https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-18864/E-GEOD-18864.raw.3.zip",
              destfile = "./data-raw/GSE18864/E-GEOD-18864.raw.3.zip") #Downloading the CEL files.
unzip("./data-raw/GSE18864/E-GEOD-18864.raw.1.zip", exdir = "./data-raw/GSE18864/GSE18864-raw/") #Extract files
unzip("./data-raw/GSE18864/E-GEOD-18864.raw.2.zip", exdir = "./data-raw/GSE18864/GSE18864-raw/") #Extract files
unzip("./data-raw/GSE18864/E-GEOD-18864.raw.3.zip", exdir = "./data-raw/GSE18864/GSE18864-raw/") #Extract files
unlink("./data-raw/GSE18864/E-GEOD-18864.raw.1.zip") #Removing the zip file
unlink("./data-raw/GSE18864/E-GEOD-18864.raw.2.zip") #Removing the zip file
unlink("./data-raw/GSE18864/E-GEOD-18864.raw.3.zip") #Removing the zip file

#Also downloading the annotations:
download.file(url = "https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-18864/E-GEOD-18864.sdrf.txt",
              destfile = "./data-raw/GSE18864/GSE18864-raw/E-GEOD-18864.sdrf.txt")


###Preparing Gene Expression Data:
##We use affy package to do RMA of the CEL files:
require(affy)
EGEOD18864 <- justRMA(celfile.path = "./data-raw/GSE18864/GSE18864-raw/")

#Extracting gene expression matrix:
EGEOD18864GEX <- exprs(EGEOD18864)

#Converting affy gene IDs to Entrez IDs:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("affy_hg_u133_plus_2","entrezgene"),
                     filters = "affy_hg_u133_plus_2", values = rownames(EGEOD18864GEX), mart = humaRt)

rownames(EGEOD18864GEX) <- ConvTabelle$entrezgene[match(rownames(EGEOD18864GEX),ConvTabelle$affy_hg_u133_plus_2)]
EGEOD18864GEX <- EGEOD18864GEX[!is.na(rownames(EGEOD18864GEX)),] #Losing 16699 probes!


###Making the Annotation:
#Annotation downloaded from ArrayExpress:
ExtraAnnot <- read.delim("./data-raw/GSE18864/GSE18864-raw/E-GEOD-18864.sdrf.txt")

#Order of annot not the same as data, correcting:
ExtraAnnot <- ExtraAnnot[match(substr(ExtraAnnot$Source.Name, start = 1, stop = 9),
                               substr(colnames(EGEOD18864GEX), start = 1, stop = 9)),]


##In Geeleher et al. paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4054092/) they only kept patients, will do the same:
EGEOD18864GEX_patient <- EGEOD18864GEX[,1:24]
ExtraAnnot_patient <- ExtraAnnot[1:24,]

#RECIST Annotation Geeleher et al. paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4054092/) used:
#These annotations were extracted from this paper ”https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4054092/#B50”'s
# code supplement (downloadable via http://genemed.uchicago.edu/~pgeeleher/cgpPrediction/paper.zip), found in the
# variable named ‘respOrd’ found in “Data/cisplatinData/cisplatinBreast.RData”. The order of annotations are corrected
# to match the order of our samples (starting with GSM467523 and going up till GSM467546).
annot <- as.character(read.table("./data-raw/GSE18864/GSE18864-raw/Geeleher_RECIST_for_E-GEOD-18864.txt",
                                 skip = 3)[,1]) # ( From Geeleher et. al.: "Clinical response are categorized as
                                                                                          # “clinical complete response” (cCR), “clinical partial response”
                                                                                          # (cPR), “stable disease” (SD) or “progressive disease” (PD) )
names(annot) <- ExtraAnnot_patient$Characteristics..miller.payne.response.


###Making and saving the Foresee object:
EGEOD18864 <- list()
class(EGEOD18864) <- "ForeseePatient"
EGEOD18864[["GeneExpression"]] <- EGEOD18864GEX_patient
EGEOD18864[["Annotation"]] <- annot
EGEOD18864[["ExtraAnnotation"]] <- ExtraAnnot_patient

devtools::use_data(EGEOD18864, overwrite = T)
