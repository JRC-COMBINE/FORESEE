# FORESEE: a Tool for the Systematic Comparison of Translational Drug Response Modeling Pipelines

The uniFied translatiOnal dRug rESponsE prEdiction platform FORESEE is an R-package that is designed to act as a scaffold in developing and benchmarking translational drug sensitivity models. The package is generally geared to utilize drug sensitivity knowledge gained in cancer cell line modeling to predict clinical therapy outcome in patients. For this purpose, FORESEE introduces a uniform data format for public cell line and patient data sets on the one hand and incorporates state-of-the-art preprocessing methods, model training algorithms and different validation techniques on the other hand. The modular implementation of these different functional elements offers the training and testing of diverse combinatorial models, which can be used to re-evaluate and improve already existing modeling pipelines, but also to develop new ones.


### Table of Contents
**[Installation](#installation)**<br>
**[License](#license)**<br>
**[Usage Instructions](#usage-instructions)**<br>
**[Test Environments](#test-environments)**<br>

### Installation

#### Full Version (Code+Data)
Because of file size restrictions on github, FORESEE on github does not contain any data. 
If you want to install the full FORESEE package, first download the package via this link: 
https://osf.io/k3pg7/download
Or you can download directly in R:
```r
download.file(url = "https://osf.io/k3pg7/download", destfile = "FORESEE_0.9.9.tar.gz")
#You can change destination folder as you wish; here, we are downloading to the current folder
```
After downloading you can install the package by:
```r
install.packages("FORESEE_0.9.9.tar.gz", repos = NULL, type = "source")
```

#### GitHub Version (Code only)
If you are only interested on the methods provided in FORESEE but not the data sets, you can install FORESEE from this repository.
Easiest way to install GitHub version of FORESEE is via <a href="https://cran.r-project.org/web/packages/devtools/">Devtools</a>.
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
The package was tested on R 3.4, on Windows 10, Mac OS X and Linux (CENTOS 7.4).
