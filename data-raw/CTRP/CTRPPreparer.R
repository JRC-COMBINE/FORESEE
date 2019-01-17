## Script that will download, process and add CTRP(v2) to this package
# by esfahani@aices.rwth-aachen.de, XX.1.2019

## Raw data were downloaded from ctd2 portal (
#https://ocg.cancer.gov/programs/ctd2/data-portal) on 1.2019 . All needed
#files are stored in a zip file, which can be downloaded and unzipped via:
download.file(url =
                "ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip",
              destfile = "data-raw/CTRP/CTRP-raw/CTRPv2.0_2015_ctd2_ExpandedDataset.zip")
unzip(zipfile = "data-raw/CTRP/CTRP-raw/CTRPv2.0_2015_ctd2_ExpandedDataset.zip",
      exdir = "data-raw/CTRP/CTRP-raw/")
unlink(c("data-raw/CTRP/CTRP-raw/CTRPv2.0_2015_ctd2_ExpandedDataset.zip", #Removing big unwanted files
         "data-raw/CTRP/CTRP-raw/v20.data.per_cpd_avg.txt",
         "data-raw/CTRP/CTRP-raw/v20.data.per_cpd_post_qc.txt",
         "data-raw/CTRP/CTRP-raw/v20.data.per_cpd_pre_qc.txt",
         "data-raw/CTRP/CTRP-raw/v20.data.per_cpd_well.txt"))
# Or you can download CTRPv2.0_2015_ctd2_ExpandedDataset.zip manually from
#CDT2 (or any other resource listed in CTRP website:
#https://portals.broadinstitute.org/ctrp/) and put its content in CTRP-raw
#folder.

## Here's a list of raw files, expected to be in data-raw/CTRP/CTRP-raw/:
#0- v20._README.txt and v20._COLUMNS.txt are not actually needed for making
#CTRP object; they are included since they contain extra information about
#other files used in making CTRP.
#1- v20.data.curves_post_qc.txt is the main CTRP file, containing the
#experiment results. From v20._README.txt, it contains 'area-under-concentration-response
#curve (AUC) sensitivity scores, curve-fit parameters, and confidence intervals
#for each cancer cell line and each compound'
#2- v20.meta.per_compound.txt: 'contextual compound information and annotation'
#(text from v20._README.txt).
#3- v20.meta.per_cell_line.txt: 'contextual cancer cell line information and annotation'
#(text from v20._README.txt).
#4- v20.meta.per_experiment.txt: 'information about experimental growth conditions,
#media, and SNP fingerprinting' (text from v20._README.txt).


###Response data curation:
##Loading up relevant files:
ResponseData <- read.table("./data-raw/CTRP/CTRP-raw/v20.data.curves_post_qc.txt",
                           header = TRUE)
CompoundInfo <- read.delim("./data-raw/CTRP/CTRP-raw/v20.meta.per_compound.txt",
                           header = TRUE)
CellInfo <- read.delim("./data-raw/CTRP/CTRP-raw/v20.meta.per_cell_line.txt",
                           header = TRUE)
ExperimentInfo <- read.delim("./data-raw/CTRP/CTRP-raw/v20.meta.per_experiment.txt",
                           header = TRUE)
##Some merging of different tables:
ExperimentInfo[["cell_name"]] <- CellInfo$ccl_name[match(ExperimentInfo$master_ccl_id,
                                                         CellInfo$master_ccl_id)]
ResponseData[["cell_name"]] <- ExperimentInfo$cell_name[match(ResponseData$experiment_id,
                                                              ExperimentInfo$experiment_id)]
ResponseData[["cpd_name"]] <- CompoundInfo$cpd_name[match(ResponseData$master_cpd_id,
                                                          CompoundInfo$master_cpd_id)]

## We rearrange drug responses (EC50, AUC, ...) into a matrix, with cell lines as rows and drugs as columns:
#Initializing drug response matrices:
EC50 <- AUC <- PPV <- matrix(NA, nrow = nlevels(ResponseData$cell_name),
                                  ncol = nlevels(ResponseData$cpd_name),
                                  dimnames = list(levels(ResponseData$cell_name),
                                                  levels(ResponseData$cpd_name)))

#We use progress package for making a more sophisticated progress bar for our loops:
require(progress)
pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                       total = nrow(ResponseData), clear = FALSE)

#These loops fill drug response matrices we initialized above:
for(KEI in seq_len(nrow(ResponseData))){
  CellLine <- ResponseData$cell_name[KEI]
  DrUg <- ResponseData$cpd_name[KEI]
  EC50[CellLine, DrUg] <- ResponseData$apparent_ec50_umol[KEI]
  AUC[CellLine, DrUg] <- ResponseData$area_under_curve[KEI]
  PPV[CellLine, DrUg] <- ResponseData$pred_pv_high_conc[KEI]
  pb$tick()
}


###Preparing 'ResponseTypes' slot:
##ResponseTypes is a data frame containing variable names and description of
#availble ForeseeCell drug response variables(Info from v20.COLUMNS.txt):
Response <- c("EC50","AUC","PPV")
Response <- cbind(Response,c("apparent_ec50_umol: apparent effective concentration at which 50% of total decline is observed",
                             "area_under_curve: integrated area under the sigmoid-fit concentration-response curve",
                             "pred_pv_high_conc: predicted percent viability at highest concentration tested"))
Response <- as.data.frame(Response);colnames(Response) <- c("Name","Description")

###Preparing 'InputTypes' slot:
##InputTypes is a data frame containing variable names (and description) of
#availble ForeseeCell input data variables (identical to CCLE):
InputTypes <- c("GeneExpression","GeneExpressionRNAseq","ProteinExpression","CNVGain","CNVLoss","Mutation")
InputTypes <- cbind(InputTypes,c("RMA-normalized gene expression measured by DNA array",
                                 "Gene expression measured by sequencing in RPKM",
                                 "Protein expression measured by Reverse Phase Protein Array",
                                 "Gaining of copy number variation, a binary matrix",
                                 "Losing of copy number variation, a binary matrix",
                                 "Mutation, a binary matrix"))
InputTypes <- as.data.frame(InputTypes);colnames(InputTypes) <- c("Name","Description")

###Preparing Meta information:
names(CellInfo)[4] <- "Site" #For being compatible with SampleSelector function in FORESEE
#Drug Info:
names(CompoundInfo)[c(2,7)] <- c("DRUG_NAME","TARGET") #For being compatible with FeatureSelector function in FORESEE

###Preparing CTRP Genetic feature data:
##Based on CTRP references, e.g. Basu et al.(Cell 2013), 'Genetic feature data
#for CCLs tested [in CTRP] can be downloaded from the Broad/Novartis CCLE portal
#(http://www.broadinstitute.org/ccle).'. Since we already prepared an object
#from CCLE data portal, we use the already prepared genetic data in CCLE object.
#We expect CCLE preparation is already done before running this script. if not,
#and/or for more info, check data-raw/CCLE/CCLEPreparer.R
# Loading up CCLE:
load("data/CCLE.rda")
# Adapting CCLE cell names to match names in CTRP data (by removing tissue name):
colnames(CCLE$GeneExpression) <- sapply(strsplit(colnames(CCLE$GeneExpression),
                                                 split = "_"), function(x) x[1])
colnames(CCLE$GeneExpressionRNAseq) <- sapply(strsplit(colnames(CCLE$GeneExpressionRNAseq),
                                                       split = "_"), function(x) x[1])
colnames(CCLE$ProteinExpression) <- sapply(strsplit(colnames(CCLE$ProteinExpression),
                                                    split = "_"), function(x) x[1])
colnames(CCLE$CNVGain) <- sapply(strsplit(colnames(CCLE$CNVGain),
                                          split = "_"), function(x) x[1])
colnames(CCLE$CNVLoss) <- sapply(strsplit(colnames(CCLE$CNVLoss),
                                          split = "_"), function(x) x[1])
colnames(CCLE$Mutation) <- sapply(strsplit(colnames(CCLE$Mutation),
                                           split = "_"), function(x) x[1])


###Making and saving the Foresee object:
CTRP <- list()
class(CTRP) <- "ForeseeCell"
CTRP[["GeneExpression"]] <- CCLE$GeneExpression
CTRP[["GeneExpressionRNAseq"]] <- CCLE$GeneExpressionRNAseq
CTRP[["ProteinExpression"]]<- CCLE$ProteinExpression
CTRP[["CNVGain"]] <- CCLE$CNVGain
CTRP[["CNVLoss"]] <- CCLE$CNVLoss
CTRP[["Mutation"]] <- CCLE$Mutation
CTRP[["EC50"]] <- EC50
CTRP[["AUC"]] <- AUC
CTRP[["PPV"]] <- PPV
CTRP[["DrugInfo"]] <- CompoundInfo
CTRP[["TissueInfo"]] <- CellInfo
CTRP[["ExperimentInfo"]] <- ExperimentInfo
CTRP[["ResponseTypes"]] <- Response
CTRP[["InputTypes"]] <- InputTypes

devtools::use_data(CTRP, overwrite = T)

