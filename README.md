# ActiveDriverWGS

ActiveDriverWGS is a driver discovery tool for analysis of whole genome sequencing data. It analyzes the mutational burden of SNVs and short INDELs in functionally defined regions of interest and retrieves regions which are significantly mutated compared to a background window.

For more information, please refer to the [ActiveDriverWGS preprint](https://www.biorxiv.org/content/early/2017/12/19/236802).

## Installation

#### devtools:
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/ActiveDriverWGSR')`

#### From source
Clone the repository: `https://github.com/reimandlab/ActiveDriverWGSR.git`
Open R in the directory you cloned the package in and run `install.packages('ActiveDriverWGSR', repos=NULL)`

## Using ActiveDriverWGS
#### NOTE: Currently ActiveDriverWGS only works on regions/mutations from hg19. Stay tuned for an update.

For a tutorial on ActiveDriverWGS, please click [here](http://htmlpreview.github.io/?http://github.com/reimandlab/ActiveDriverWGSR/blob/master/vignettes/ActiveDriverWGSR.html).

[test](http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html)
