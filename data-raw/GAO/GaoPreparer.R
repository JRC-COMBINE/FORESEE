## Script that will (download,) process and add GAO Xenograft dataset to this package
## data published in https://www.ncbi.nlm.nih.gov/pubmed/26479923
# by esfahani@aices.rwth-aachen.de, 15.3.2018

## We downloaded GAO Xenograft dataset as supplement files of the paper
#"High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response"
#by Gao et. al. published in 2015 in nature medicine accessible via the link https://www.nature.com/articles/nm.3954
# Paper itself is locked behind a pay wall, but the supplement we used is freely available and can be downloaded
#with the line:
download.file(url = "https://media.nature.com/original/nature-assets/nm/journal/v21/n11/extref/nm.3954-S2.xlsx",
              destfile = "./data-raw/GAO/GAO-raw/nm.3954-S2.xlsx")

###Preparing Gene Expression Data based on RNA-seq:
#Importing and cleaning up:
require(openxlsx)
RNAseq_fpkm <- read.xlsx(xlsxFile = "./data-raw/GAO/GAO-raw/nm.3954-S2.xlsx", sheet = 1)
RNAseq_fpkmMat <- as.matrix(RNAseq_fpkm[,-1])

##We move the fpkm values into logarithmic (base 2) scale so that we'd be having semi-normaly distributed values, necessary for linear models
#We cut all lower that 1 values, and replace them with 1.
RNAseq_fpkmMat[RNAseq_fpkmMat < 1] <- 1
RNAseq_fpkmMat <- log2(RNAseq_fpkmMat)
rownames(RNAseq_fpkmMat) <- RNAseq_fpkm$Sample

#Feature identifiers are in gene symbols, we convert them to Entrez IDs:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol", values = rownames(RNAseq_fpkmMat), mart = humaRt)
rownames(RNAseq_fpkmMat) <- ConvTabelle$entrezgene[match(rownames(RNAseq_fpkmMat),ConvTabelle$hgnc_symbol)]
RNAseq_fpkmMat <- RNAseq_fpkmMat[!is.na(rownames(RNAseq_fpkmMat)),]


###Preparing Copy number data(SNP array):
#Importing and cleaning up:
copy_number <- read.xlsx(xlsxFile = "./data-raw/GAO/GAO-raw/nm.3954-S2.xlsx", sheet = 2)
copy_numberMat <- as.matrix(copy_number[,-1])
rownames(copy_numberMat) <- copy_number$Sample

#Removing 'ArmLevelCNScore' and 'FocalCNScore' (I didn't use them in anyway, e.g. normalize other values based on them):
copy_numberMat <- copy_numberMat[-c(1,2),]

#Converting gene symbols to Entrez IDs:
ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol", values = rownames(copy_numberMat), mart = humaRt)
rownames(copy_numberMat) <- ConvTabelle$entrezgene[match(rownames(copy_numberMat),ConvTabelle$hgnc_symbol)]
copy_numberMat <- copy_numberMat[!is.na(rownames(copy_numberMat)),]


###Preparing Copy number and Mutation data(Targeted sequencing by SureSelect Human All Exon V4, Agilent):
#Importing and cleaning up:
mut_and_cn2 <- read.xlsx(xlsxFile = "./data-raw/GAO/GAO-raw/nm.3954-S2.xlsx", sheet = 3)
#Based on unique(mut_and_cn2$Category) we have two different kind of amplification and 3 types of mutation,
#but we aim to make one binary matrix for amplification, one for deletion and one for mutation, so we assume
#all different types of mutation and amplifications are the same:
CNVGainIndx <- mut_and_cn2$Category=="Amp5" | mut_and_cn2$Category=="Amp8"
CNVLossIndx <- mut_and_cn2$Category=="Del0.8"
MutationIndx <- mut_and_cn2$Category=="MutKnownFunctional" | mut_and_cn2$Category=="MutLikelyFunctional" | mut_and_cn2$Category=="MutNovel"

CNVGainTmp <- mut_and_cn2[CNVGainIndx,]
CNVLossTmp <- mut_and_cn2[CNVLossIndx,]
MutationTmp <- mut_and_cn2[MutationIndx,]

#After separating mut_and_cn2 into CNV gain, CNV loss and mutation data.frames, we arrage
#all three into three binary matices with genes as rows and samples as columns:

#Starting with CNVGain matrix:
CNVGain <- matrix(data = 0, nrow = length(unique(CNVGainTmp$Entrez)), ncol = length(unique(CNVGainTmp$Sample)),
                  dimnames = list(unique(CNVGainTmp$Entrez),unique(CNVGainTmp$Sample)))
#We use progress package for making a more sophisticated progress bar for our loops:
require(progress)
pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                       total = nrow(CNVGainTmp), clear = FALSE)
# Filling the binary matrix:
for(K in 1:nrow(CNVGainTmp)){
  CNVGain[CNVGainTmp$Entrez[K],CNVGainTmp$Sample[K]] <- 1
  pb$tick()
}
# There were duplications in CNVGainTmp (for example check mut_and_cn2[751:752,]), bug in data?? ignore for now.

#Continuing with CNVLoss matrix:
CNVLoss <- matrix(data = 0, nrow = length(unique(CNVLossTmp$Entrez)), ncol = length(unique(CNVLossTmp$Sample)),
                  dimnames = list(unique(CNVLossTmp$Entrez),unique(CNVLossTmp$Sample)))
#We use progress package for making a more sophisticated progress bar for our loops:
require(progress)
pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                       total = nrow(CNVLossTmp), clear = FALSE)
# Filling the binary matrix:
for(K in 1:nrow(CNVLossTmp)){
  CNVLoss[CNVLossTmp$Entrez[K],CNVLossTmp$Sample[K]] <- 1
  pb$tick()
}
# Again there were duplications in data, bug in data?? ignore for now.

#Finishing with Mutation matrix:
Mutation <- matrix(data = 0, nrow = length(unique(MutationTmp$Entrez)), ncol = length(unique(MutationTmp$Sample)),
                  dimnames = list(unique(MutationTmp$Entrez),unique(MutationTmp$Sample)))
#We use progress package for making a more sophisticated progress bar for our loops:
require(progress)
pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                       total = nrow(MutationTmp), clear = FALSE)
# Filling the binary matrix:
for(K in 1:nrow(MutationTmp)){
  Mutation[MutationTmp$Entrez[K],MutationTmp$Sample[K]] <- 1
  pb$tick()
}
# Again there were duplications in data, bug in data?? ignore for now.


###Response data curation:
#Importing and cleaning up:
PCT_metrics <- read.xlsx(xlsxFile = "./data-raw/GAO/GAO-raw/nm.3954-S2.xlsx", sheet = 5)
#There are single and double drug treatment in this dataset; we separate them into two different metrics:
PCT_metrics_single <- PCT_metrics[PCT_metrics$Treatment.type=="single",]
PCT_metrics_combo <- PCT_metrics[PCT_metrics$Treatment.type=="combo",]

##We make a matrix for each and every response measure, with model in rows and treatment in columns:
#First for single treatment:
PCT_metrics_single_BestResponse <- PCT_metrics_single_Day_BestResponse <- PCT_metrics_single_BestAvgResponse <-
  PCT_metrics_single_Day_BestAvgResponse <- PCT_metrics_single_TimeToDouble <- PCT_metrics_single_Day_Last <-
  PCT_metrics_single_ResponseCategory <-
  matrix(NA,nrow = length(unique(PCT_metrics_single$Model)), ncol = length(unique(PCT_metrics_single$Treatment)),
         dimnames = list(unique(PCT_metrics_single$Model),unique(PCT_metrics_single$Treatment)))
for(K in 1:nrow(PCT_metrics_single)){
  PCT_metrics_single_BestResponse[PCT_metrics_single$Model[K],PCT_metrics_single$Treatment[K]] <- PCT_metrics_single$BestResponse[K]
  PCT_metrics_single_Day_BestResponse[PCT_metrics_single$Model[K],PCT_metrics_single$Treatment[K]] <- PCT_metrics_single$Day_BestResponse[K]
  PCT_metrics_single_BestAvgResponse[PCT_metrics_single$Model[K],PCT_metrics_single$Treatment[K]] <- PCT_metrics_single$BestAvgResponse[K]
  PCT_metrics_single_Day_BestAvgResponse[PCT_metrics_single$Model[K],PCT_metrics_single$Treatment[K]] <- PCT_metrics_single$Day_BestAvgResponse[K]
  PCT_metrics_single_TimeToDouble[PCT_metrics_single$Model[K],PCT_metrics_single$Treatment[K]] <- PCT_metrics_single$TimeToDouble[K]
  PCT_metrics_single_Day_Last[PCT_metrics_single$Model[K],PCT_metrics_single$Treatment[K]] <- PCT_metrics_single$Day_Last[K]
  PCT_metrics_single_ResponseCategory[PCT_metrics_single$Model[K],PCT_metrics_single$Treatment[K]] <- PCT_metrics_single$ResponseCategory[K]
}

#Then for double treatments:
PCT_metrics_combo_BestResponse <- PCT_metrics_combo_Day_BestResponse <- PCT_metrics_combo_BestAvgResponse <-
  PCT_metrics_combo_Day_BestAvgResponse <- PCT_metrics_combo_TimeToDouble <- PCT_metrics_combo_Day_Last <-
  PCT_metrics_combo_ResponseCategory <-
  matrix(NA,nrow = length(unique(PCT_metrics_combo$Model)), ncol = length(unique(PCT_metrics_combo$Treatment)),
         dimnames = list(unique(PCT_metrics_combo$Model),unique(PCT_metrics_combo$Treatment)))
for(K in 1:nrow(PCT_metrics_combo)){
  PCT_metrics_combo_BestResponse[PCT_metrics_combo$Model[K],PCT_metrics_combo$Treatment[K]] <- PCT_metrics_combo$BestResponse[K]
  PCT_metrics_combo_Day_BestResponse[PCT_metrics_combo$Model[K],PCT_metrics_combo$Treatment[K]] <- PCT_metrics_combo$Day_BestResponse[K]
  PCT_metrics_combo_BestAvgResponse[PCT_metrics_combo$Model[K],PCT_metrics_combo$Treatment[K]] <- PCT_metrics_combo$BestAvgResponse[K]
  PCT_metrics_combo_Day_BestAvgResponse[PCT_metrics_combo$Model[K],PCT_metrics_combo$Treatment[K]] <- PCT_metrics_combo$Day_BestAvgResponse[K]
  PCT_metrics_combo_TimeToDouble[PCT_metrics_combo$Model[K],PCT_metrics_combo$Treatment[K]] <- PCT_metrics_combo$TimeToDouble[K]
  PCT_metrics_combo_Day_Last[PCT_metrics_combo$Model[K],PCT_metrics_combo$Treatment[K]] <- PCT_metrics_combo$Day_Last[K]
  PCT_metrics_combo_ResponseCategory[PCT_metrics_combo$Model[K],PCT_metrics_combo$Treatment[K]] <- PCT_metrics_combo$ResponseCategory[K]
}


###Preparing Meta information:
#Importing Meta Info to as tissue information:
PCT_raw <- read.xlsx(xlsxFile = "./data-raw/GAO/GAO-raw/nm.3954-S2.xlsx", sheet = 4)
names(PCT_raw)[2] <- "Site" #For being compatible with SampleSelector function in FORESEE

#Drug Info:
DInfo <- PCT_metrics[,2:4]
names(DInfo)[1:2] <- c("DRUG_NAME","TARGET") #For being compatible with FeatureSelector function in FORESEE

###Preparing 'ResponseTypes' slot:
##ResponseTypes is a data frame containing variable names and description of availble ForeseeCell drug response variables:
Response <- c("BestResponse","Day_BestResponse","BestAvgResponse","Day_BestAvgResponse","TimeToDouble","Day_Last","ResponseCategory",
              "BestResponseCombo","Day_BestResponseCombo","BestAvgResponseCombo","Day_BestAvgResponseCombo","TimeToDoubleCombo","Day_LastCombo","ResponseCategoryCombo")
Response <- cbind(Response,c("The BestResponse was the minimum value of ΔVolt(tumor volume change vs. time t to its baseline) for t ≥ 10 d (For single treatment)",
                             "The day of the the BestResponse (For single treatment)",
                             "The best average of ΔVolt(tumor volume change vs. time t to its baseline) from t = 0 to t (For single treatment)",
                             "The day of the the best average of ΔVolt (For single treatment)",
                             "Time till tumor volume doubled (For single treatment)",
                             "Duration of treatment in days (For single treatment)",
                             "The criteria for response (mRECIST) were adapted from RECIST criteria and defined as follows (applied in this order):
                             mCR, BestResponse < −95% and BestAvgResponse < −40%; mPR, BestResponse < −50% and BestAvgResponse < −20%;
                             mSD, BestResponse < 35% and BestAvgResponse < 30%; mPD, not otherwise categorized. (For single treatment)",
                             "The BestResponse was the minimum value of ΔVolt(tumor volume change vs. time t to its baseline) for t ≥ 10 d (For double treatment)",
                             "The day of the the BestResponse (For double treatment)",
                             "The best average of ΔVolt(tumor volume change vs. time t to its baseline) from t = 0 to t (For double treatment)",
                             "The day of the the best average of ΔVolt (For double treatment)",
                             "Time till tumor volume doubled (For double treatment)",
                             "Duration of treatment in days (For double treatment)",
                             "The criteria for response (mRECIST) were adapted from RECIST criteria and defined as follows (applied in this order):
                             mCR, BestResponse < −95% and BestAvgResponse < −40%; mPR, BestResponse < −50% and BestAvgResponse < −20%;
                             mSD, BestResponse < 35% and BestAvgResponse < 30%; mPD, not otherwise categorized. (For double treatment)"))
Response <- as.data.frame(Response);colnames(Response) <- c("Name","Description")

###Preparing 'InputTypes' slot:
##InputTypes is a data frame containing variable names (and description) of availble ForeseeCell input data variables:
InputTypes <- c("GeneExpression","SNP6","CNVGain","CNVLoss","Mutation")
InputTypes <- cbind(InputTypes,c("Log2 value of gene expression measured by sequencing in fpkm",
                                 "Copy-number variation measured by SNP array (Affymetrix genome-wide human SNP Array 6.0 chip)",
                                 "Gaining of copy number variation, a binary matrix, measured by Targeted sequencing (SureSelect Human All Exon V4, Agilent)",
                                 "Losing of copy number variation, a binary matrix, measured by Targeted sequencing (SureSelect Human All Exon V4, Agilent)",
                                 "Mutation, a binary matrix, measured by Targeted sequencing (SureSelect Human All Exon V4, Agilent)"))
InputTypes <- as.data.frame(InputTypes);colnames(InputTypes) <- c("Name","Description")



###Making and saving the Foresee object:
GAO <- list()
class(GAO) <- "ForeseeCell" ##We assign Xenografts as Cell objects!
GAO[["GeneExpression"]] <- RNAseq_fpkmMat
GAO[["SNP6"]] <- copy_numberMat
GAO[["CNVGain"]] <- CNVGain
GAO[["CNVLoss"]] <- CNVLoss
GAO[["Mutation"]] <- Mutation

GAO[["BestResponse"]] <- PCT_metrics_single_BestResponse
GAO[["Day_BestResponse"]] <- PCT_metrics_single_Day_BestResponse
GAO[["BestAvgResponse"]] <- PCT_metrics_single_BestAvgResponse
GAO[["Day_BestAvgResponse"]] <- PCT_metrics_single_Day_BestAvgResponse
GAO[["TimeToDouble"]] <- PCT_metrics_single_TimeToDouble
GAO[["Day_Last"]] <- PCT_metrics_single_Day_Last
GAO[["ResponseCategory"]] <- PCT_metrics_single_ResponseCategory

GAO[["BestResponseCombo"]] <- PCT_metrics_combo_BestResponse
GAO[["Day_BestResponseCombo"]] <- PCT_metrics_combo_Day_BestResponse
GAO[["BestAvgResponseCombo"]] <- PCT_metrics_combo_BestAvgResponse
GAO[["Day_BestAvgResponseCombo"]] <- PCT_metrics_combo_Day_BestAvgResponse
GAO[["TimeToDoubleCombo"]] <- PCT_metrics_combo_TimeToDouble
GAO[["Day_LastCombo"]] <- PCT_metrics_combo_Day_Last
GAO[["ResponseCategoryCombo"]] <- PCT_metrics_combo_ResponseCategory

GAO[["DrugInfo"]] <- DInfo
GAO[["TissueInfo"]] <- PCT_raw

GAO[["ResponseTypes"]] <- Response
GAO[["InputTypes"]] <- InputTypes

devtools::use_data(GAO, overwrite = T)

