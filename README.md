# ActiveDriverWGS

ActiveDriverWGS is a cancer driver discovery tool for analysis of somatic mutations derived from whole genome sequencing. It works on protein-coding sequences as well as various non-coding sequences (non-coding RNAs, promoters, enhancers, to name a few). ActiveDriverWGS is based on a statistical model that analyzes the mutational burden of SNVs and short indels in functionally defined elements of interest. It then retrieves elements which are significantly mutated compared to a background sequence window, its nucleotide sequence context and trinucleotide mutation signatures.

Two major kinds of input are required by ActiveDriverWGS: a set of genomic regions of interest defined as a BED12 file, and a set of somatic mutations in a cohort of tumor samples in a tab-separated text file.

For more information, please refer to the [publication](https://www.cell.com/molecular-cell/fulltext/S1097-2765(19)30957-8).

Helen Zhu*, Liis Uuskula-Reimand*, Keren Isaev*, Lina Wadi, Azad Alizada, Shimin Shuai, Vincent Huang, Dike Aduluso-Nwaobasi, Marta Paczkowska, Diala Abd-Rabbo, Oliver Ocsenas, Minggao Liang, J. Drew Thompson, Yao Li, Luyao Ruan, Michal Krassowski, Irakli Dzneladze, Jared T. Simpson, Mathieu Lupien, Lincoln D. Stein, Paul C. Boutros, Michael D. Wilson, Jüri Reimand. <b>Candidate Cancer Driver Mutations in Distal Regulatory Elements and Long-Range Chromatin Interaction Networks</b>. <i>Molecular Cell</i> (2020), https://doi.org/10.1016/j.molcel.2019.12.027

## Installation

#### devtools:
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/ActiveDriverWGSR')`

#### From source
Clone the repository: `https://github.com/reimandlab/ActiveDriverWGSR.git`
Open R in the directory you cloned the package in and run `install.packages('ActiveDriverWGS', repos=NULL)`

## Using ActiveDriverWGS
#### NOTE: Currently ActiveDriverWGS only works on regions/mutations from hg19. Stay tuned for an update.

For a tutorial on ActiveDriverWGS, please click [here](http://htmlpreview.github.io/?https://github.com/reimandlab/ActiveDriverWGSR/blob/master/vignettes/ActiveDriverWGSR.html).


