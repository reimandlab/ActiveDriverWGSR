# ActiveDriverWGS

<b>Note (Sept 3rd 2022): ActiveDriverWGS (1.2.0) has been released. It includes an improved trinucleotide quantification in genomic elements and mutations. Testing for depletion of mutations in genomic elements has been enabled. Improved estimates of expected mutations have been added. See NEWS.md for detailed updates.</b>

ActiveDriverWGS is a cancer driver discovery tool for analysis of somatic mutations derived from whole genome sequencing. It works on protein-coding sequences as well as various non-coding sequences (non-coding RNAs, promoters, enhancers, etc). ActiveDriverWGS is based on a statistical model that analyzes the mutational burden of SNVs and short indels in functionally defined elements of interest. It then retrieves elements that are significantly mutated compared to their background sequence windows, nucleotide sequence contexts and trinucleotide mutation signatures. ActiveDriverWGS is implemented as an R package and is available from GitHub and CRAN.

Two major kinds of input are required by ActiveDriverWGS: a set of genomic regions of interest defined as a BED12 file, and a set of somatic mutations in a cohort of tumor samples in a tab-separated text file. The default reference genome is human (hg19)and the updated version of ActiveDriverWGS also supports the human genome (hg38) and the mouse genome (mm9, mm10). 



## Installation

#### From the CRAN repository in R:
`install.packages('ActiveDriverWGS')`

#### Using devtools in R:
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/ActiveDriverWGSR')`

#### From source:
Clone the repository: `https://github.com/reimandlab/ActiveDriverWGSR.git`
Open R in the directory you cloned the package in and run `install.packages('ActiveDriverWGS', repos = NULL)`

## Using ActiveDriverWGS

The basic use case of ActiveDriverWGS involves a set of genomic elements (each including one or more sub-elements such as exons) and a set of somatic mutations in a dataset of whole cancer genomes. Single-nucleotide variants (SNVs) and insertions-deletions (indels) are analysed. Patient IDs are required since the analysis considers at most one mutation per element per patient for driver detection. 

```R

library(ActiveDriverWGS)

##
# load test data from the package and print a few lines of each dataset
##

data("cll_mutations")
head(cll_mutations, n = 2)
#          chr      pos1      pos2 ref alt       patient
# 779701  chr6  96651182  96651183  AA   T 001-0002-03TD
# 779702 chr10 106556005 106556005   C   T 001-0002-03TD

data("cancer_genes")
head(cancer_genes, n = 2)
#      chr   start     end       id
# 648 chr1 2488103 2488172 TNFRSF14
# 649 chr1 2489164 2489273 TNFRSF14

##
# quick example: test a few genes
##
some_genes = c("ATM", "MYD88", "NOTCH1")
elements_to_test = cancer_genes[cancer_genes$id %in% some_genes,]

results = ActiveDriverWGS(mutations = cll_mutations, elements = elements_to_test)

results
#       id  pp_element element_muts_obs element_muts_exp element_enriched pp_site
# 1    ATM 0.002092488                2       0.11059736             TRUE       1
# 2  MYD88 0.003309630                1       0.01337728             TRUE       1
# 3 NOTCH1 0.004964277                2       0.16562752             TRUE       1
#   site_muts_obs site_muts_exp site_enriched fdr_element fdr_site
# 1            NA            NA            NA 0.004964277        1
# 2            NA            NA            NA 0.004964277        1
# 3            NA            NA            NA 0.004964277        1
#   has_site_mutations
# 1
# 2
# 3

##
# load coordinates of elements from a BED12 file (a function for BED4 format is also available). 
##

elements = prepare_elements_from_BED12(
		system.file("extdata", "chr17.coding_regions.bed", 
				package = "ActiveDriverWGS"))

head(elements, n = 2)
#     chr start   end                                             id
# 1 chr17  6006  6168 gc19_pc.cds::gencode::DOC2B::ENSG00000272636.1
# 2 chr17 11205 11332 gc19_pc.cds::gencode::DOC2B::ENSG00000272636.1

results = ActiveDriverWGS(mutations = cll_mutations, elements = elements)

```

## Detailed tutorial

For a tutorial on ActiveDriverWGS, please click [here](http://htmlpreview.github.io/?https://github.com/reimandlab/ActiveDriverWGSR/blob/master/doc/ActiveDriverWGSR.html).

## Publication

For more information, please refer to the [publication](https://doi.org/10.1016/j.molcel.2019.12.027).

Helen Zhu*, Liis Uuskula-Reimand*, Keren Isaev*, Lina Wadi, Azad Alizada, Shimin Shuai, Vincent Huang, Dike Aduluso-Nwaobasi, Marta Paczkowska, Diala Abd-Rabbo, Oliver Ocsenas, Minggao Liang, J. Drew Thompson, Yao Li, Luyao Ruan, Michal Krassowski, Irakli Dzneladze, Jared T. Simpson, Mathieu Lupien, Lincoln D. Stein, Paul C. Boutros, Michael D. Wilson, Jüri Reimand. <b>Candidate Cancer Driver Mutations in Distal Regulatory Elements and Long-Range Chromatin Interaction Networks</b>. <i>Molecular Cell</i> (2020), https://doi.org/10.1016/j.molcel.2019.12.027



