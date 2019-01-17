#' Genomics of Drug Sensitivity in Cancer or GDSC
#'
#' Genomics of Drug Sensitivity in Cancer, or GDSC for short, is one of the ForeseeCell
#' datasets available in the FORESEE package. All files related to the GDSC
#' dataset were downloaded on 25.4.2018 from \url{https://www.cancerrxgene.org/downloads}.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseeCell object is very similar to the list data type in R programming
#' language; it is a data structure which includes different data types, and can be
#' indexed using double brackets or dollar sign, for example, ForeseeCell\$variable1 or
#' ForeseeCell[["variable1"]] or ForeseeCell[[1]].
#'
#' Available components in GDSC are:
#'
#' \describe{
#'   \item{GeneExpression}{RMA normalized DNA array values were converted into an R
#'   matrix, with genes in rows and cell lines in columns. Column names that were
#'   originally COSMIC IDs were converted to cell line names, and row names were
#'   converted from Ensembl gene IDs to Entrez IDs using biomaRt.}
#'   \item{CNVGain, CNVLoss, Mutation and Methylation}{Four different binary matrices,
#'   with genes in rows and cell lines in columns, all extracted from supplement
#'   Table S3B of the paper Iorio et. al. 2016. Gene names were converted from
#'   symbols to Entrez IDs.}
#'   \item{LN_IC50, AUC ,RMSE ,Z_SCORE ,MAX_CONC_MICROMOLAR , MIN_CONC_MICROMOLAR}{
#'   Response data were all in one data frame. We rearranged drug responses into
#'   different matrices, with cell lines as rows and drugs as columns.}
#'   \item{DrugInfo}{A data frame with extra information about the included drugs
#'   in the dataset. Columns 'DRUG_NAME' and 'TARGET' from this component are used
#'   in FeatureSelector.ontology() and FeatureSelector.pathway() and are needed if
#'   the user wants to use any of the mentioned FeatureSelector methods, where the
#'   pipeline uses only the gene names for training the model that are contained in
#'   the ontology or pathway associated with the chosen drug.}
#'   \item{TissueInfo}{A data frame with tissue-related information about cell lines
#'   in the dataset. Column 'Site' of this component is used to extract relevant
#'   samples that user set by assigning a value to 'TrainingTissue' input in
#'   ForeseeTrain(). Hence, this component is needed if user wants to use a
#'   specific tissue for 'TrainingTissue'.}
#'   \item{CelllineInfo}{A data frame with extra information about cell lines
#'   in the dataset.}
#'   \item{InputTypes}{InputTypes is a data frame with two columns of 'Name'
#'   and 'Description', which provide the names of all components in the object
#'   that can be used as input data (in ForeseeTrain for example), and description
#'   for each input data.}
#'   \item{ResponseTypes}{ResponseTypes is another two-column data frame with a
#'   'Name' column, providing the names of all components in the object that are
#'   a measure of drug activity and can be used as response variable (called
#'   'CellResponseType' in ForeseeTrain) and a 'Description' column for each response variable.}
#' }
#' @source \url{https://www.cancerrxgene.org/downloads}
"GDSC"

#' Broad Institute Cancer Cell Line Encyclopedia or CCLE
#'
#' Broad Institute Cancer Cell Line Encyclopedia, or CCLE for short, is a cell line
#' dataset included as a ForeseeCell instance. All relevant
#' files for the CCLE object were downloaded on May 2018 from
#' \url{https://portals.broadinstitute.org/ccle/data}.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseeCell object is very similar to the list data type in R programming
#' language; it is a data structure which includes different data types, and can be
#' indexed using double brackets or dollar sign, for example, ForeseeCell\$variable1 or
#' ForeseeCell[["variable1"]] or ForeseeCell[[1]].
#'
#' Available components in CCLE are:
#'
#' \describe{
#'   \item{GeneExpression}{A matrix, with gene Entrez IDs in rows and cell lines in
#'   columns, containing RMA-normalized gene expression profiles measured by DNA array.}
#'   \item{GeneExpressionRNAseq}{A matrix, with genes identified by Entrez IDs in rows
#'   and cell lines in columns, containing gene expression profiles measured by RNA-seq,
#'   in Reads Per Kilobase of transcript per Million (RPKM).}
#'   \item{Mutation, CNVGain and CNVLoss}{Three different binary matrices, with genes
#'   in rows (names converted to Entrez IDs) and cell lines in columns, pointing toward
#'   mutation, gaining copy number variation and losing copy number variation respectively.}
#'   \item{ProteinExpression}{Another matrix, with genes in rows and cell lines in columns,
#'   containing values of protein abundance measured by Reverse Phase Protein Array.
#'   Gene corresponding to measured protein (Target genes of antibody) in each row is
#'   identified by its Entrez ID. Since there are duplications in measured genes, we
#'   used the mean value of the duplicated genes as the final value.}
#'   \item{IC50, EC50, ActArea and Amax}{From CCLE response data, four matrices
#'   were build with cell lines in rows and drugs in columns. Descriptions about
#'   each measured value are available in ResponseTypes component.
#'   \item{DrugInfo}{A data frame with extra information about the included drugs
#'   in the dataset. Columns 'DRUG_NAME' and 'TARGET' from this component are used
#'   in FeatureSelector.ontology() and FeatureSelector.pathway() and are needed if
#'   the user wants to use any of the mentioned FeatureSelector methods, where the
#'   pipeline uses only the gene names for training the model that are contained in
#'   the ontology or pathway associated with the chosen drug.}
#'   \item{TissueInfo}{A data frame with tissue-related information about cell lines
#'   in the dataset. Column 'Site' of this component is used to extract relevant
#'   samples that user set by assigning a value to 'TrainingTissue' input in
#'   ForeseeTrain(). Hence, this component is needed if user wants to use a
#'   specific tissue for 'TrainingTissue'.}
#'   \item{InputTypes}{InputTypes is a data frame with two columns of 'Name'
#'   and 'Description', which provide the names of all components in the object
#'   that can be used as input data (in ForeseeTrain for example), and description
#'   for each input data.}
#'   \item{ResponseTypes}{ResponseTypes is another two-column data frame with a
#'   'Name' column, providing the names of all components in the object that are
#'   a measure of drug activity and can be used as response variable (called
#'   'CellResponseType' in ForeseeTrain) and a 'Description' column for each response variable.}
#' }
#' @source \url{https://portals.broadinstitute.org/ccle/data}
"CCLE"

#' Cancer Therapeutics Response Portal or CTRP (v2)
#'
#' Cancer Therapeutics Response Portal, or CTRP for short, is a cell line
#' dataset included as a ForeseeCell instance. All relevant
#' files for the CTRP object were downloaded on Jan. 2019 from
#' \url{https://ocg.cancer.gov/programs/ctd2/data-portal}.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseeCell object is very similar to the list data type in R
#' programming language; it is a data structure which includes different data
#' types, and can be indexed using double brackets or dollar sign, for example,
#' ForeseeCell\$variable1 or ForeseeCell[["variable1"]] or ForeseeCell[[1]].
#'
#' Available components in CTRP are:
#'
#' \describe{
#'   \item{GeneExpression}{A matrix, with gene Entrez IDs in rows and cell lines in
#'   columns, containing RMA-normalized gene expression profiles measured by DNA array.}
#'   \item{GeneExpressionRNAseq}{A matrix, with genes identified by Entrez IDs in rows
#'   and cell lines in columns, containing gene expression profiles measured by RNA-seq,
#'   in Reads Per Kilobase of transcript per Million (RPKM).}
#'   \item{Mutation, CNVGain and CNVLoss}{Three different binary matrices, with genes
#'   in rows (names converted to Entrez IDs) and cell lines in columns, pointing toward
#'   mutation, gaining copy number variation and losing copy number variation respectively.}
#'   \item{ProteinExpression}{Another matrix, with genes in rows and cell lines in columns,
#'   containing values of protein abundance measured by Reverse Phase Protein Array.
#'   Gene corresponding to measured protein (Target genes of antibody) in each row is
#'   identified by its Entrez ID. Since there are duplications in measured genes, we
#'   used the mean value of the duplicated genes as the final value.}
#'   \item{EC50, AUC and PPV}{From CTRP response data, three matrices
#'   were build with cell lines in rows and drugs in columns. Descriptions about
#'   each measured value are available in ResponseTypes component.
#'   \item{DrugInfo}{A data frame with extra information about the tested
#'   'small-molecule's of CTRP. Columns 'DRUG_NAME' and 'TARGET' from this component are used
#'   in FeatureSelector.ontology() and FeatureSelector.pathway() and are needed if
#'   the user wants to use any of the mentioned FeatureSelector methods, where the
#'   pipeline uses only the gene names for training the model that are contained in
#'   the ontology or pathway associated with the chosen drug.}
#'   \item{TissueInfo}{A data frame with tissue-related information about cell lines
#'   in the dataset. Column 'Site' of this component is used to extract relevant
#'   samples that user set by assigning a value to 'TrainingTissue' input in
#'   ForeseeTrain(). Hence, this component is needed if user wants to use a
#'   specific tissue for 'TrainingTissue'.}
#'   \item{InputTypes}{InputTypes is a data frame with two columns of 'Name'
#'   and 'Description', which provide the names of all components in the object
#'   that can be used as input data (in ForeseeTrain for example), and description
#'   for each input data.}
#'   \item{ResponseTypes}{ResponseTypes is another two-column data frame with a
#'   'Name' column, providing the names of all components in the object that are
#'   a measure of drug activity and can be used as response variable (called
#'   'CellResponseType' in ForeseeTrain) and a 'Description' column for each response variable.}
#' }
#' @source \url{https://portals.broadinstitute.org/ctrp.v2.1}
"CTRP"

#' DAEMEN Breast Cancer Cell Line Data Set
#'
#' DAEMEN ForeseeCell contains the data used in Daemen et. al. 2013 publication. We
#' downloaded the data used for modeling, as provided by the paper with the link
#' \url{https://www.synapse.org/#!Synapse:syn2179898}.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseeCell object is very similar to the list data type in R programming
#' language; it is a data structure which includes different data types, and can be
#' indexed using double brackets or dollar sign, for example, ForeseeCell\$variable1 or
#' ForeseeCell[["variable1"]] or ForeseeCell[[1]].
#'
#' Available components in DAEMEN are:
#'
#' \describe{
#'   \item{GeneExpression}{A matrix, with gene Entrez IDs in rows and cell lines in
#'   columns, containing RMA-normalized gene expression profiles measured by DNA array.}
#'   \item{GeneExpressionRNAseq}{Matrix of counts based on RNA-seq technology. We
#'   transformed the seq values into logarithmic (base 2) scale for having
#'   semi-normal distributed values, which is necessary for linear modeling.
#'   As a prerequisite for transforming to log-scale, we replaced all values
#'   lower than 1 with 1.}
#'   \item{Methylation and SNP6}{Directly imported into matrices from downloaded files.}
#'   \item{GI50}{A matrix with cell lines in rows and drugs in
#'   columns. A description about GI50 can be found in ResponseTypes component.}
#'   \item{TissueInfo}{A data frame with tissue-related information about cell lines
#'   in the dataset. Column 'Site' of this component is used to extract relevant
#'   samples that user set by assigning a value to 'TrainingTissue' input in
#'   ForeseeTrain(). Hence, this component is needed if user wants to use a
#'   specific tissue for 'TrainingTissue'.}
#'   \item{InputTypes}{InputTypes is a data frame with two columns of 'Name'
#'   and 'Description', which provide the names of all components in the object
#'   that can be used as input data (in ForeseeTrain for example), and description
#'   for each input data.}
#'   \item{ResponseTypes}{ResponseTypes is another two-column data frame with a
#'   'Name' column, providing the names of all components in the object that are
#'   a measure of drug activity and can be used as response variable (called
#'   'CellResponseType' in ForeseeTrain) and a 'Description' column for each response variable.}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/24176112}
"DAEMEN"

#' GAO Xenograft Data Set
#'
#' GAO is one of the two xenograft data sets included in FORESEE. The data is downloaded as supplement files of
#' Gao et. al. 2015 (https://www.ncbi.nlm.nih.gov/pubmed/26479923), which are freely available and can be downloaded
#' via https://media.nature.com/original/nature-assets/nm/journal/v21/n11/extref/nm.3954-S2.xlsx.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseeCell object is very similar to the list data type in R programming
#' language; it is a data structure which includes different data types, and can be
#' indexed using double brackets or dollar sign, for example, ForeseeCell\$variable1 or
#' ForeseeCell[["variable1"]] or ForeseeCell[[1]].
#'
#' Available components in GAO are:
#'
#' \describe{
#'   \item{GeneExpression}{A matrix of RNA-seq values in FPKM (Fragments Per Kilobase of transcript per Million).
#'   The matrix contains genes as rows and samples as columns, with Entrez IDs in rownames. FPKM values were
#'   transformed into logarithmic scale (base 2) for having semi-normal distributed values, which is necessary
#'   for linear modeling.}
#'   \item{SNP6}{Copy number data measured by SNP array (Affymetrix genome-wide human SNP Array 6.0 chip)}
#'   \item{Mutation, CNVGain and CNVLoss}{Binary matrices, pointing toward mutations, gaining copy number
#'   variations and losing copy number variations respectively.}
#'   \item{BestResponse, BestResponseCombo, TimeToDouble, ...}{This data set includes 14 different response
#'   matrices, all of which have samples in rows and drugs in columns. List and description of these
#'   response matrices can be found in ResponseTypes component.}
#'   \item{TissueInfo}{A data frame with tissue-related information about cell lines
#'   in the dataset. Column 'Site' of this component is used to extract relevant
#'   samples that user set by assigning a value to 'TrainingTissue' input in
#'   ForeseeTrain(). Hence, this component is needed if user wants to use a
#'   specific tissue for 'TrainingTissue'.}
#'   \item{DrugInfo}{A data frame with extra information about the included drugs
#'   in the dataset. Columns 'DRUG_NAME' and 'TARGET' from this component are used
#'   in FeatureSelector.ontology() and FeatureSelector.pathway() and are needed if
#'   the user wants to use any of the mentioned FeatureSelector methods, where the
#'   pipeline uses only the gene names for training the model that are contained in
#'   the ontology or pathway associated with the chosen drug.}
#'   \item{InputTypes}{InputTypes is a data frame with two columns of 'Name'
#'   and 'Description', which provide the names of all components in the object
#'   that can be used as input data (in ForeseeTrain for example), and description
#'   for each input data.}
#'   \item{ResponseTypes}{ResponseTypes is another two-column data frame with a
#'   'Name' column, providing the names of all components in the object that are
#'   a measure of drug activity and can be used as response variable (called
#'   'CellResponseType' in ForeseeTrain) and a 'Description' column for each response variable.}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/26479923}
"GAO"

#' WITKIEWICZ Xenograft Data Set
#'
#' WITKIEWICZ is the other xenograft data set included in FORESEE. This data set is from a study by Witkiewcz et al. 2016,
#' studying Pancreatic ductal adenocarcinoma (PDAC) drug response.
#' Data used in building the WITKIEWICZ ForeseeCell are two excel files, which are included as supplements in the original
#' paper and the GEO data set GSE84023, which includes the RNA-seq gene expression relevant to the paper.
#'
#' REMEMBER listDrugs(OBJ = WITKIEWICZ) or listInputOptions(FunArgument = "DrugName", OBJ = WITKIEWICZ)) only
#' return drugs of single treatment (acceptable as DrugName in ForeseeTrain when CellResponseType="AUC"), if
#' you want to have CellResponseType="AUCCombo" in ForeseeTrain, list of possible "DrugName" values can be
#' listed by colnames(WITKIEWICZ$AUCCombo).
#'
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#'
#'
#' @format A ForeseeCell object is very similar to the list data type in R programming
#' language; it is a data structure which includes different data types, and can be
#' indexed using double brackets or dollar sign, for example, ForeseeCell\$variable1 or
#' ForeseeCell[["variable1"]] or ForeseeCell[[1]].
#'
#' Available components in WITKIEWICZ are:
#'
#' \describe{
#'   \item{GeneExpression}{Matrix of gene expressions measured by RNA-seq. We used the already processed
#'   data available on GEO. Based on GSE84023 page on GEO, this is the processing pipeline they used:
#'   "Illumina Casava1.7 software used for basecalling. Sequenced reads were trimmed for adaptor sequence,
#'   mapped to hg19 genome using bowtie TopHat. Counts per gene was obtained using HTseq counts and normalized
#'   using edgeR package in R. Genome_build: hg19. files_format_andy_content: tab-delimited text file include
#'   matrix of normalized log counts per million for each sample."
#'   We averaged over all samples from the same patient.}
#'   \item{AUC and AUCCombo}{Response data were imported from supplement excel files, and then formatted as a
#'   matrix with samples in rows and drugs in columns. More information about these two response matrices can be
#'   found in ResponseTypes component.
#'   In AUC, samples from the same patient are averaged.}
#'   \item{DrugInfo}{A data frame with extra information about the included drugs
#'   in the dataset. Columns 'DRUG_NAME' and 'TARGET' from this component are used
#'   in FeatureSelector.ontology() and FeatureSelector.pathway() and are needed if
#'   the user wants to use any of the mentioned FeatureSelector methods, where the
#'   pipeline uses only the gene names for training the model that are contained in
#'   the ontology or pathway associated with the chosen drug.}
#'   \item{InputTypes}{InputTypes is a data frame with two columns of 'Name'
#'   and 'Description', which provide the names of all components in the object
#'   that can be used as input data (in ForeseeTrain for example), and description
#'   for each input data.}
#'   \item{ResponseTypes}{ResponseTypes is another two-column data frame with a
#'   'Name' column, providing the names of all components in the object that are
#'   a measure of drug activity and can be used as response variable (called
#'   'CellResponseType' in ForeseeTrain) and a 'Description' column for each response variable.}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/27498862}
"WITKIEWICZ"

#' GSE6434 Patient Data Set
#'
#' GSE6434 is a gene expression response dataset of 24 patients to Docetaxel treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in GSE6434 are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were downloaded from GEO, and normalized using RMA from affy package.}
#'   \item{Annotation}{is a logical or numeric vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation), e.g. when Annotation is a logical vector, names(Annotation) provides information
#'   about what True and False in Annotation mean in terms of patient response.}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6434}
"GSE6434"

#' EGEOD18864 Patient Data Set
#'
#' EGEOD18864 (identical to GSE18864 on GEO) is a gene expression response dataset of 24 patients to Cisplatin treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in EGEOD18864 are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were downloaded from Array Express, and normalized using RMA from affy package.}
#'   \item{Annotation}{is a logical or numeric vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation), e.g. when Annotation is a logical vector, names(Annotation) provides information
#'   about what True and False in Annotation mean in terms of patient response.}
#'   \item{ExtraAnnotation}{is a data frame including all annotations that was contained in the original
#'   patient data set. This component is not used in the FORESEE pipeline, but is included
#'   for the user (e.g. to divide a patient data set into sub groups based on ExtraAnnotation for better modeling).}
#' }
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-18864/}
"EGEOD18864"

#' GSE33072_erlotinib Patient Data Set
#'
#' GSE33072_erlotinib (subset of GSE33072 on GEO) is a gene expression response dataset of 25 patients to Erlotinib treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in GSE33072_erlotinib are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were downloaded from GEO, and normalized using RMA from affy package.}
#'   \item{Annotation}{is a logical or numeric vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation), e.g. when Annotation is a logical vector, names(Annotation) provides information
#'   about what True and False in Annotation mean in terms of patient response.}
#'   \item{ExtraAnnotation}{is a data frame including all annotations that was contained in the original
#'   patient data set. This component is not used in the FORESEE pipeline, but is included
#'   for the user (e.g. to divide a patient data set into sub groups based on ExtraAnnotation for better modeling).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33072}
"GSE33072_erlotinib"

#' GSE33072_sorafenib Patient Data Set
#'
#' GSE33072_sorafenib (subset of GSE33072 on GEO) is a gene expression response dataset of 39 patients to Sorafenib treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in GSE33072_sorafenib are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were downloaded from GEO, and normalized using RMA from affy package.}
#'   \item{Annotation}{is a logical or numeric vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation), e.g. when Annotation is a logical vector, names(Annotation) provides information
#'   about what True and False in Annotation mean in terms of patient response.}
#'   \item{ExtraAnnotation}{is a data frame including all annotations that was contained in the original
#'   patient data set. This component is not used in the FORESEE pipeline, but is included
#'   for the user (e.g. to divide a patient data set into sub groups based on ExtraAnnotation for better modeling).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33072}
"GSE33072_sorafenib"

#' GSE9782_GPL96_bortezomib Patient Data Set
#'
#' GSE9782_GPL96_bortezomib (subset of GSE9782 on GEO, measured by Affymetrix Human Genome U133A Array) is a gene expression
#' response dataset of 169 patients to Bortezomib treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in GSE9782_GPL96_bortezomib are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were not available at GEO, so we downloaded the MAS5.0 normalized values from GEO and
#'   transformed the data to a logarithmic scale (base 2).}
#'   \item{Annotation}{is a logical or numeric vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation), e.g. when Annotation is a logical vector, names(Annotation) provides information
#'   about what True and False in Annotation mean in terms of patient response.}
#'   \item{ExtraAnnotation}{is a data frame including all annotations that was contained in the original
#'   patient data set. This component is not used in the FORESEE pipeline, but is included
#'   for the user (e.g. to divide a patient data set into sub groups based on ExtraAnnotation for better modeling).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9782}
"GSE9782_GPL96_bortezomib"

#' GSE9782_GPL96_dexamethasone Patient Data Set
#'
#' GSE9782_GPL96_dexamethasone (subset of GSE9782 on GEO, measured by Affymetrix Human Genome U133A Array) is a gene expression
#' response dataset of 70 patients to Dexamethasone treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in GSE9782_GPL96_dexamethasone are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were not available at GEO, so we downloaded the MAS5.0 normalized values from GEO and
#'   transformed the data to a logarithmic scale (base 2).}
#'   \item{Annotation}{is a logical or numeric vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation), e.g. when Annotation is a logical vector, names(Annotation) provides information
#'   about what True and False in Annotation mean in terms of patient response.}
#'   \item{ExtraAnnotation}{is a data frame including all annotations that was contained in the original
#'   patient data set. This component is not used in the FORESEE pipeline, but is included
#'   for the user (e.g. to divide a patient data set into sub groups based on ExtraAnnotation for better modeling).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9782}
"GSE9782_GPL96_dexamethasone"

#' GSE9782_GPL97_bortezomib Patient Data Set
#'
#' GSE9782_GPL97_bortezomib (subset of GSE9782 on GEO, measured by Affymetrix Human Genome U133B Array) is a gene expression
#' response dataset of 169 patients to Bortezomib treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in GSE9782_GPL97_bortezomib are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were not available at GEO, so we downloaded the MAS5.0 normalized values from GEO and
#'   transformed the data to a logarithmic scale (base 2).}
#'   \item{Annotation}{is a logical or numeric vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation), e.g. when Annotation is a logical vector, names(Annotation) provides information
#'   about what True and False in Annotation mean in terms of patient response.}
#'   \item{ExtraAnnotation}{is a data frame including all annotations that was contained in the original
#'   patient data set. This component is not used in the FORESEE pipeline, but is included
#'   for the user (e.g. to divide a patient data set into sub groups based on ExtraAnnotation for better modeling).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9782}
"GSE9782_GPL97_bortezomib"

#' GSE9782_GPL97_dexamethasone Patient Data Set
#'
#' GSE9782_GPL97_dexamethasone (subset of GSE9782 on GEO, measured by Affymetrix Human Genome U133B Array) is a gene expression
#' response dataset of 70 patients to Dexamethasone treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in GSE9782_GPL97_dexamethasone are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were not available at GEO, so we downloaded the MAS5.0 normalized values from GEO and
#'   transformed the data to a logarithmic scale (base 2).}
#'   \item{Annotation}{is a logical or numeric vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation), e.g. when Annotation is a logical vector, names(Annotation) provides information
#'   about what True and False in Annotation mean in terms of patient response.}
#'   \item{ExtraAnnotation}{is a data frame including all annotations that was contained in the original
#'   patient data set. This component is not used in the FORESEE pipeline, but is included
#'   for the user (e.g. to divide a patient data set into sub groups based on ExtraAnnotation for better modeling).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9782}
"GSE9782_GPL97_dexamethasone"

#' GSE51373 Patient Data Set
#'
#' GSE51373 is a gene expression response dataset of 28 serous epithelial ovarian
#' cancer patients to combinatory platin-based chemo plus Paclitaxel treatment.
#' You can check the data vignette for more information (browseVignettes(package = "FORESEE")).
#'
#' @format A ForeseePatient is a data structure having components of different data types,
#' that can be indexed using double brackets or dollar sign (for example
#' ForeseePatient\$variable1 or ForeseePatient[["variable1"]] or ForeseePatient[[1]]),
#' similar to a list data type in R programming language.
#'
#' Available components in GSE51373 are:
#'
#' \describe{
#'   \item{GeneExpression}{is a matrix, with genes in rows and patients in columns. Entrez IDs are
#'   saved in 'rownames' of the matrix and patient identifiers in 'colnames'.
#'
#'   Raw CEL files were downloaded from GEO, and normalized using RMA from affy package.}
#'   \item{Annotation}{is a logical vector indicating the patient response to a drug. Extra information
#'   is provided in names(Annotation).}
#'   \item{ExtraAnnotation}{is a data frame including all annotations that was contained in the original
#'   patient data set. This component is not used in the FORESEE pipeline, but is included
#'   for the user (e.g. to divide a patient data set into sub groups based on ExtraAnnotation for better modeling).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51373}
"GSE51373"

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


