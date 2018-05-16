## Script that will download, process and add GSE9782 to FORESEE package
# by esfahani@aices.rwth-aachen.de, 23.2.2018

##Raw CEL files of GSE9782 are not available in Gene Expression Omnibus (GEO), check the
#link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9782.
##Since only normalized data was available, we had to go with normalized data (Which could
#be problematic, specially because the normalization is MAS5 and not RMA, but since
#Geeleher et. al. in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4054092/ used this
#dataset and we aim to replicate their work, we used this dataset too.)
## You have to downlaod Series Matrix Files (in txt format) and put them in GSE9782-raw,
#You can do this process automatically by running the next lines:
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9782/matrix//GSE9782-GPL96_series_matrix.txt.gz",
              destfile = "./data-raw/GSE9782/GSE9782-raw/GSE9782-GPL96_series_matrix.txt.gz")
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9782/matrix//GSE9782-GPL97_series_matrix.txt.gz",
              destfile = "./data-raw/GSE9782/GSE9782-raw/GSE9782-GPL97_series_matrix.txt.gz")


##You need the GEOquery package (https://bioconductor.org/packages/release/bioc/html/GEOquery.html) to import GSE files:
library(GEOquery)

###Preparing Gene Expression Data:
##Importing GSE files:
GSE9782_GPL96 <- getGEO(filename = "./data-raw/GSE9782/GSE9782-raw//GSE9782-GPL96_series_matrix.txt.gz", getGPL = F)
GSE9782_GPL97 <- getGEO(filename = "./data-raw/GSE9782/GSE9782-raw//GSE9782-GPL97_series_matrix.txt.gz", getGPL = F)

#Extracting gene expression matrices:
GSE9782_GPL96Mat <- exprs(GSE9782_GPL96)
GSE9782_GPL97Mat <- exprs(GSE9782_GPL97)

#Converting affy gene IDs to Entrez IDs:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("affy_hg_u133a","entrezgene"),
                     filters = "affy_hg_u133a", values = rownames(GSE9782_GPL96Mat), mart = humaRt)

rownames(GSE9782_GPL96Mat) <- ConvTabelle$entrezgene[match(rownames(GSE9782_GPL96Mat),ConvTabelle$affy_hg_u133a)]
GSE9782_GPL96Mat <- GSE9782_GPL96Mat[!is.na(rownames(GSE9782_GPL96Mat)),]

ConvTabelle <- getBM(attributes = c("affy_hg_u133b","entrezgene"),
                     filters = "affy_hg_u133b", values = rownames(GSE9782_GPL97Mat), mart = humaRt)

rownames(GSE9782_GPL97Mat) <- ConvTabelle$entrezgene[match(rownames(GSE9782_GPL97Mat),ConvTabelle$affy_hg_u133b)]
GSE9782_GPL97Mat <- GSE9782_GPL97Mat[!is.na(rownames(GSE9782_GPL97Mat)),]


###Dividing GSE9782 into sub-data-sets:
##In GSE9782 there are two sub-datasets measured by U133A and U133B Arrays (GPL97 and GPL97 respectively), each of these
#sub-datasets can be also divided into more sub-sets based on the treatment of the patients.
##We already divided GSE9782 into GSE9782_GPL96 and GSE9782_GPL97, we will also divide each onto two subtypes of
# 'responders and non-responders to treatment with bortezomib' and 'responders and non-responders to treatment with dexamethasone'.
#(Some samples will be discarded with this sub-setting, we used these categorization we wanted to replicate Geeleher et. al.
#in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4054092/).

##Starting with bortezomib samples in GPL96:
#Extracting sub-gene-epression-matrix:
GSE9782_GPL96Mat_bortezomib_NR_R <- GSE9782_GPL96Mat[,GSE9782_GPL96@phenoData@data$characteristics_ch1.1 == "treatment = PS341" &
                                                       GSE9782_GPL96@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
#Making the Annotation:
GSE9782_GPL96Mat_bortezomib_NR_R_ExtraAnnot <- GSE9782_GPL96@phenoData@data[GSE9782_GPL96@phenoData@data$characteristics_ch1.1 == "treatment = PS341" &
                                                                  GSE9782_GPL96@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE",]
GSE9782_GPL96Mat_bortezomib_Annotation_Labels <- GSE9782_GPL96@phenoData@data$characteristics_ch1.7[GSE9782_GPL96@phenoData@data$characteristics_ch1.1 == "treatment = PS341" &
                                                                                               GSE9782_GPL96@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
GSE9782_GPL96Mat_bortezomib_Annotation <- GSE9782_GPL96@phenoData@data$characteristics_ch1.8[GSE9782_GPL96@phenoData@data$characteristics_ch1.1 == "treatment = PS341" &
                                                                                                      GSE9782_GPL96@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
#Since we have responder/nonresponder, we make a 'logical vector' as annotation:
GSE9782_GPL96Mat_bortezomib_Annotation <- GSE9782_GPL96Mat_bortezomib_Annotation == "PGx_Responder = R"
names(GSE9782_GPL96Mat_bortezomib_Annotation) <- GSE9782_GPL96Mat_bortezomib_Annotation_Labels


##Continuing with bortezomib samples in GPL97:
#Extracting sub-gene-epression-matrix:
GSE9782_GPL97Mat_bortezomib_NR_R <- GSE9782_GPL97Mat[,GSE9782_GPL97@phenoData@data$characteristics_ch1.1 == "treatment = PS341" &
                                                       GSE9782_GPL97@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
#Making the Annotation:
GSE9782_GPL97Mat_bortezomib_NR_R_ExtraAnnot <- GSE9782_GPL97@phenoData@data[GSE9782_GPL97@phenoData@data$characteristics_ch1.1 == "treatment = PS341" &
                                                                              GSE9782_GPL97@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE",]
GSE9782_GPL97Mat_bortezomib_Annotation_Labels <- GSE9782_GPL97@phenoData@data$characteristics_ch1.7[GSE9782_GPL97@phenoData@data$characteristics_ch1.1 == "treatment = PS341" &
                                                                                                      GSE9782_GPL97@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
GSE9782_GPL97Mat_bortezomib_Annotation <- GSE9782_GPL97@phenoData@data$characteristics_ch1.8[GSE9782_GPL97@phenoData@data$characteristics_ch1.1 == "treatment = PS341" &
                                                                                               GSE9782_GPL97@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
#Since we have responder/nonresponder, we make a 'logical vector' as annotation:
GSE9782_GPL97Mat_bortezomib_Annotation <- GSE9782_GPL97Mat_bortezomib_Annotation == "PGx_Responder = R"
names(GSE9782_GPL97Mat_bortezomib_Annotation) <- GSE9782_GPL97Mat_bortezomib_Annotation_Labels


##Continuing with dexamethasone samples in GPL96:
#Extracting sub-gene-epression-matrix:
GSE9782_GPL96Mat_dexamethasone_NR_R <- GSE9782_GPL96Mat[,GSE9782_GPL96@phenoData@data$characteristics_ch1.1 == "treatment = Dex" &
                                                       GSE9782_GPL96@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
#Making the Annotation:
GSE9782_GPL96Mat_dexamethasone_NR_R_ExtraAnnot <- GSE9782_GPL96@phenoData@data[GSE9782_GPL96@phenoData@data$characteristics_ch1.1 == "treatment = Dex" &
                                                                              GSE9782_GPL96@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE",]
GSE9782_GPL96Mat_dexamethasone_Annotation_Labels <- GSE9782_GPL96@phenoData@data$characteristics_ch1.7[GSE9782_GPL96@phenoData@data$characteristics_ch1.1 == "treatment = Dex" &
                                                                                                      GSE9782_GPL96@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
GSE9782_GPL96Mat_dexamethasone_Annotation <- GSE9782_GPL96@phenoData@data$characteristics_ch1.8[GSE9782_GPL96@phenoData@data$characteristics_ch1.1 == "treatment = Dex" &
                                                                                               GSE9782_GPL96@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
#Since we have responder/nonresponder, we make a 'logical vector' as annotation:
GSE9782_GPL96Mat_dexamethasone_Annotation <- GSE9782_GPL96Mat_dexamethasone_Annotation == "PGx_Responder = R"
names(GSE9782_GPL96Mat_dexamethasone_Annotation) <- GSE9782_GPL96Mat_dexamethasone_Annotation_Labels


##Finishing with dexamethasone samples in GPL97:
#Extracting sub-gene-epression-matrix:
GSE9782_GPL97Mat_dexamethasone_NR_R <- GSE9782_GPL97Mat[,GSE9782_GPL97@phenoData@data$characteristics_ch1.1 == "treatment = Dex" &
                                                       GSE9782_GPL97@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
#Making the Annotation:
GSE9782_GPL97Mat_dexamethasone_NR_R_ExtraAnnot <- GSE9782_GPL97@phenoData@data[GSE9782_GPL97@phenoData@data$characteristics_ch1.1 == "treatment = Dex" &
                                                                                 GSE9782_GPL97@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE",]
GSE9782_GPL97Mat_dexamethasone_Annotation_Labels <- GSE9782_GPL97@phenoData@data$characteristics_ch1.7[GSE9782_GPL97@phenoData@data$characteristics_ch1.1 == "treatment = Dex" &
                                                                                                      GSE9782_GPL97@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
GSE9782_GPL97Mat_dexamethasone_Annotation <- GSE9782_GPL97@phenoData@data$characteristics_ch1.8[GSE9782_GPL97@phenoData@data$characteristics_ch1.1 == "treatment = Dex" &
                                                                                               GSE9782_GPL97@phenoData@data$characteristics_ch1.8 != "PGx_Responder = IE"]
#Since we have responder/nonresponder, we make a 'logical vector' as annotation:
GSE9782_GPL97Mat_dexamethasone_Annotation <- GSE9782_GPL97Mat_dexamethasone_Annotation == "PGx_Responder = R"
names(GSE9782_GPL97Mat_dexamethasone_Annotation) <- GSE9782_GPL97Mat_dexamethasone_Annotation_Labels



###Unlike Geeleher et. al. in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4054092/, we moved
#values of GSE9782 into log-scale because in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9782
#they mentioned "... Data was normalized to a Ttimmed mean of 15o and is NOT log transformed. ...". All
#other datasets in FORESEE are RMA normalized so they are in log-scale, so we do the same for GSE9782:
GSE9782_GPL96Mat_bortezomib_NR_R_loged <- log2(GSE9782_GPL96Mat_bortezomib_NR_R)
GSE9782_GPL97Mat_bortezomib_NR_R_loged <- log2(GSE9782_GPL97Mat_bortezomib_NR_R)
GSE9782_GPL96Mat_dexamethasone_NR_R_loged <- log2(GSE9782_GPL96Mat_dexamethasone_NR_R)
GSE9782_GPL97Mat_dexamethasone_NR_R_loged <- log2(GSE9782_GPL97Mat_dexamethasone_NR_R)



###Making and saving the Foresee object(s):
GSE9782_GPL96_bortezomib <- list()
class(GSE9782_GPL96_bortezomib) <- "ForeseePatient"
GSE9782_GPL96_bortezomib[["GeneExpression"]] <- GSE9782_GPL96Mat_bortezomib_NR_R_loged
GSE9782_GPL96_bortezomib[["Annotation"]] <- GSE9782_GPL96Mat_bortezomib_Annotation
GSE9782_GPL96_bortezomib[["ExtraAnnotation"]] <- GSE9782_GPL96Mat_bortezomib_NR_R_ExtraAnnot

devtools::use_data(GSE9782_GPL96_bortezomib, overwrite = T)

GSE9782_GPL97_bortezomib <- list()
class(GSE9782_GPL97_bortezomib) <- "ForeseePatient"
GSE9782_GPL97_bortezomib[["GeneExpression"]] <- GSE9782_GPL97Mat_bortezomib_NR_R_loged
GSE9782_GPL97_bortezomib[["Annotation"]] <- GSE9782_GPL97Mat_bortezomib_Annotation
GSE9782_GPL97_bortezomib[["ExtraAnnotation"]] <- GSE9782_GPL97Mat_bortezomib_NR_R_ExtraAnnot

devtools::use_data(GSE9782_GPL97_bortezomib, overwrite = T)

GSE9782_GPL96_dexamethasone <- list()
class(GSE9782_GPL96_dexamethasone) <- "ForeseePatient"
GSE9782_GPL96_dexamethasone[["GeneExpression"]] <- GSE9782_GPL96Mat_dexamethasone_NR_R_loged
GSE9782_GPL96_dexamethasone[["Annotation"]] <- GSE9782_GPL96Mat_dexamethasone_Annotation
GSE9782_GPL96_dexamethasone[["ExtraAnnotation"]] <- GSE9782_GPL96Mat_dexamethasone_NR_R_ExtraAnnot

devtools::use_data(GSE9782_GPL96_dexamethasone, overwrite = T)

GSE9782_GPL97_dexamethasone <- list()
class(GSE9782_GPL97_dexamethasone) <- "ForeseePatient"
GSE9782_GPL97_dexamethasone[["GeneExpression"]] <- GSE9782_GPL97Mat_dexamethasone_NR_R_loged
GSE9782_GPL97_dexamethasone[["Annotation"]] <- GSE9782_GPL97Mat_dexamethasone_Annotation
GSE9782_GPL97_dexamethasone[["ExtraAnnotation"]] <- GSE9782_GPL97Mat_dexamethasone_NR_R_ExtraAnnot

devtools::use_data(GSE9782_GPL97_dexamethasone, overwrite = T)



