## Notice: This repository is deprecated!
Due to file size restrictions of GitHub, we moved FORESEE to [a GitLab repository, hosted by RWTH Aachen university](https://git.rwth-aachen.de/jrc-combine/foresee), and continued support and development there. 

Therefore, for accessing the last version of FORESEE, use the GitLab version here: https://git.rwth-aachen.de/jrc-combine/foresee.

# FORESEE: a Tool for the Systematic Comparison of Translational Drug Response Modeling Pipelines

The uniFied translatiOnal dRug rESponsE prEdiction platform FORESEE is an R-package that is designed to act as a scaffold in developing and benchmarking translational drug sensitivity models. The package is generally geared to utilize drug sensitivity knowledge gained in cancer cell line modeling to predict clinical therapy outcome in patients. For this purpose, FORESEE introduces a uniform data format for public cell line and patient data sets on the one hand and incorporates state-of-the-art preprocessing methods, model training algorithms and different validation techniques on the other hand. The modular implementation of these different functional elements offers the training and testing of diverse combinatorial models, which can be used to re-evaluate and improve already existing modeling pipelines, but also to develop new ones.


### Table of Contents
**[Installation](#installation)**<br>
**[License](#license)**<br>
**[Usage Instructions](#usage-instructions)**<br>
**[Test Environments](#test-environments)**<br>

### Installation

#### Full Version (Code+Data)
FORESEE is available on GitHub and OSF (https://osf.io/rf6qk/ or DOI 10.17605/OSF.IO/RF6QK). Because of file size restrictions on github, FORESEE on github does not contain any data. 
If you want to install the full FORESEE package, you have to install the OSF version. At the time of writing this README, the last available version of FORESEE on OSF is 1.1.1. You can install this version by firstly downloading the source package via this link:
https://osf.io/h43uy/download
or directly downloading it in R:
```r
download.file(url = "https://osf.io/h43uy/download", destfile = "FORESEE_1.1.1.tar.gz")
#You can change destination folder as you wish; here, we are downloading to the current folder
```
After downloading you can install the source package by:
```r
install.packages("FORESEE_1.1.1.tar.gz", repos = NULL, type = "source")
```

#### GitHub Version (Code only)
If you are only interested on the methods provided in FORESEE but not the data sets, you can install FORESEE from this repository.
The easiest way to install the GitHub version of FORESEE is via <a href="https://cran.r-project.org/web/packages/devtools/">Devtools</a>.
After installing Devtools from cran, you can install FORESEE by:
```r
devtools::install_github(repo = "JRC-COMBINE/FORESEE", build_vignettes = TRUE)
```

### License
Distributed under GNU General Public License v3.0. See the accompanying [license](https://github.com/JRC-COMBINE/FORESEE/blob/master/LICENSE) file or the copy at https://www.gnu.org/licenses/gpl-3.0.html.

### Usage Instructions

The package contains a vignette to explain the preparation of data contained in the package accessible via
```r
vignette(topic="DataOverview")
```
and a vignette to explain some of the functional routines that can be implemented using FORESEE accessible via 
```r
vignette(topic="introduction-to-FORESEE")
```

You can additionally access the vignettes mentioned above in your browser by:
```r
browseVignettes(package = "FORESEE")
```

### Test Environments
The package was tested with R 3.4 on Windows 10, Mac OS X and Linux (CentOS 7.4).
