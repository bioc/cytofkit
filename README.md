cytofkit: an integrated mass cytometry data analysis pipeline
============

### cytofkit

This package is designed to facilitate the analysis workflow of mass cytometry data with automatic subset identification and mapping of cellular progression. Both command line and a GUI client are provided for executing the workflow easily.

### Installation

To install the package from Bioconductor, use:

``` r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("cytofkit")
```

To install the latest version from the github repository, use:

``` r
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("JinmiaoChenLab/cytofkit")
```

### Usage

- [cytofkit: Analysis Pipeline](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_workflow.html)    
- [cytofkit: Quick Start](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_example.html)   
- [cytofkit: ShinyAPP Tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_shinyAPP.html)    
- [shiny APP](https://chenhao.shinyapps.io/cytofkitShinyAPP/) link:  https://chenhao.shinyapps.io/cytofkitShinyAPP/