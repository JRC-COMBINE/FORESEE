# FORESEE: a Tool for the Systematic Comparison of Translational Drug Response Modeling Pipelines


uniFied translatiOnal dRug rESponsE prEdiction platform FORESEE is a versatile open source software, implemented as R-package, that is designed to act as a scaffold in developing and benchmarking such computational drug response models. FORESEE not only introduces a uniform data format for public cell line and patient data sets, but also establishes a standardized environment for drug response prediction pipelines, incorporating state of the art preprocessing methods, model training algorithms and different validation techniques. A modular implementation of the different elements of the modeling pipeline facilitates a straightforward development of different combinatorial models, which can be used to re-evaluate and improve already existing modeling pipelines, as well as to develop and benchmark new ones.

We recommended to install <a href="https://github.com/JRC-COMBINE/FORESEEData">FORESEEData</a> package alongside this package.

#### Table of Contents
**[Installation Instructions](#installation-instructions)**<br>
**[Usage Instructions](#usage-instructions)**<br>
**[Troubleshooting](#troubleshooting)**<br>
**[Compatibility](#compatibility)**<br>
**[Notes and Miscellaneous](#notes-and-miscellaneous)**<br>
**[Next Steps, Credits, Feedback, License](#next-steps)**<br>

### Installation Instructions
#### Installing via Devtools (Recommended method):
Easiest way to install FORESEE is via <a href="https://cran.r-project.org/web/packages/devtools/">Devtools</a>.
After installing Devtools from cran, you can install FORESEE by:
```r
devtools::install_github(repo = "JRC-COMBINE/FORESEE")
```

#### Alternative installation methods (Manual download):
...
```Shell
cd [Your desired directory]
git clone https://github.com/JRC-COMBINE/FORESEE.git
```

```Shell
R -e "devtools::install_local('./FORESEE/')"
```

### Usage Instructions
Check the vignette of the package for instruction and tutorials:

```r
browseVignettes(package = "FORESEE")
```

### Troubleshooting
### Compatibility
### Notes and Miscellaneous
### Building the Extension Bundles
### Next Steps, Credits, Feedback, License
