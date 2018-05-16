## Script that will download, process and add GSE33072 to FORESEE package
# by esfahani@aices.rwth-aachen.de, 27.2.2018


##Raw CEL files of GSE33072 must be downloaded from Gene Expression Omnibus (GEO) via the
#link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33072, and untared in GSE33072-raw folder.
##You can do this process automatically by running the next lines:
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE33072&format=file",
              destfile = "./data-raw/GSE33072/GSE33072_RAW.tar") #Downloading the CEL files.
untar("./data-raw/GSE33072/GSE33072_RAW.tar", exdir = "./data-raw/GSE33072/GSE33072-raw/") #Extract files
unlink("./data-raw/GSE33072/GSE33072_RAW.tar") #Removing the tar file


###Preparing Gene Expression Data:
##We use affy package to do RMA of the CEL files:
require(affy)
GSE33072 <- justRMA(celfile.path = "./data-raw/GSE33072/GSE33072-raw/")

#Extracting gene expression matrix:
GSE33072GEX <- exprs(GSE33072)

#Converting affy gene IDs to Entrez IDs:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("affy_hugene_1_0_st_v1","entrezgene"),
                     filters = "affy_hugene_1_0_st_v1", values = rownames(GSE33072GEX), mart = humaRt)

rownames(GSE33072GEX) <- ConvTabelle$entrezgene[match(rownames(GSE33072GEX),ConvTabelle$affy_hugene_1_0_st_v1)]
GSE33072GEX <- GSE33072GEX[!is.na(rownames(GSE33072GEX)),]

###Making the Annotation:
##Annotation on GEO is a mess, we curated the annotation on GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33072)
#and we use our own curated annotaion:
ExtraAnnot <- readRDS("./data-raw/GSE33072/GSE33072-raw/GSE33072_patientdata_sorted.RDS")

##We divide GSE33072 into sub-datasets, based on the treatment of the patients:
#there are two treatment groups: erlotinib and sorafenib (low sample number groups are neglected; check table(ExtraAnnot$treatment))
#erlotinib:
erlotinibIndx <- which(ExtraAnnot$treatment == "erlotinib")
GSE33072GEX_erlotinib <- GSE33072GEX[,erlotinibIndx]
annot_erlotinib <- ExtraAnnot$pfsm_progression.free_survival_time_months[erlotinibIndx]
names(annot_erlotinib) <- ExtraAnnot$pfsc_progression.free.survival.status[erlotinibIndx]
ExtraAnnot_erlotinib <- ExtraAnnot[erlotinibIndx,]

#sorafenib:
sorafenibIndx <- which(ExtraAnnot$treatment == "sorafenib")
GSE33072GEX_sorafenib <- GSE33072GEX[,sorafenibIndx]
annot_sorafenib <- ExtraAnnot$pfsm_progression.free_survival_time_months[sorafenibIndx]
names(annot_sorafenib) <- ExtraAnnot$pfsc_progression.free.survival.status[sorafenibIndx]
ExtraAnnot_sorafenib <- ExtraAnnot[sorafenibIndx,]



###Making and saving the Foresee object(s):
GSE33072_erlotinib <- list()
class(GSE33072_erlotinib) <- "ForeseePatient"
GSE33072_erlotinib[["GeneExpression"]] <- GSE33072GEX_erlotinib
GSE33072_erlotinib[["Annotation"]] <- annot_erlotinib
GSE33072_erlotinib[["ExtraAnnotation"]] <- ExtraAnnot_erlotinib

devtools::use_data(GSE33072_erlotinib, overwrite = T)


GSE33072_sorafenib <- list()
class(GSE33072_sorafenib) <- "ForeseePatient"
GSE33072_sorafenib[["GeneExpression"]] <- GSE33072GEX_sorafenib
GSE33072_sorafenib[["Annotation"]] <- annot_sorafenib
GSE33072_sorafenib[["ExtraAnnotation"]] <- ExtraAnnot_sorafenib

devtools::use_data(GSE33072_sorafenib, overwrite = T)

