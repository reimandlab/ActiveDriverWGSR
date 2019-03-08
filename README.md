# ActiveDriverWGS

ActiveDriverWGS is a cancer driver discovery tool for analysis of somatic mutations derived from whole genome sequencing. It works on protein-coding sequences as well as various non-coding sequences (non-coding RNAs, promoters, enhancers, to name a few). ActiveDriverWGS is based on a statistical model that analyzes the mutational burden of SNVs and short indels in functionally defined elements of interest. It then retrieves elements which are significantly mutated compared to a background sequence window, its nucleotide sequence context and trinucleotide mutation signatures.

Two major kinds of input are required by ActiveDriverWGS: a set of genomic regions of interest defined as a BED12 file, and a set of somatic mutations in a cohort of tumor samples in a tab-separated text file.

For more information, please refer to the [ActiveDriverWGS preprint](https://www.biorxiv.org/content/early/2017/12/19/236802).

## Installation

#### devtools:
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/ActiveDriverWGS')`

#### From source
Clone the repository: `https://github.com/reimandlab/ActiveDriverWGS.git`
Open R in the directory you cloned the package in and run `install.packages('ActiveDriverWGS', repos=NULL)`

## Using ActiveDriverWGS
#### NOTE: Currently ActiveDriverWGS only works on regions/mutations from hg19. Stay tuned for an update.

For a tutorial on ActiveDriverWGS, please click [here](http://htmlpreview.github.io/?http://github.com/reimandlab/ActiveDriverWGS/blob/master/inst/doc/ActiveDriverWGS.html).
