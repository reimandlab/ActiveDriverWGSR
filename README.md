# ActiveDriverWGSR

ActiveDriverWGS is a driver discovery tool for analysis of whole genome sequencing data. It analyzes the mutational burden of SNVs and short INDELs in functionally defined regions of interest and retrieves regions which are significantly mutated compared to a background window.

For more information, please refer to the ActiveDriverWGS publication. https://www.biorxiv.org/content/early/2017/12/19/236802

## Installation

#### devtools:
Using the R package `devtools`, run
`devtools::install_github('https://github.com/reimandlab/ActiveDriverWGSR')`

#### From source
Clone the repository: `https://github.com/reimandlab/ActiveDriverWGSR.git`
Open R in the directory you cloned the package in and run `install.packages('ActiveDriverWGSR', repos=NULL)`

## Using ActiveDriverWGSR
#### NOTE: Currently ActiveDriverWGS only works on regions/mutations from GRCh37. Stay tuned for an update.
ActiveDriverWGS requires three input files:

#### 1. Mutations
6 Column Tab Separated File (Chromosome, Start, End, Reference, Alternate, Patient ID)

```
chr17	7577058	7577058	C	A	patient_1
chr17	7577094	7577094	G	A	patient_2
chr17	7578239	7578239	C	A	patient_3
chr17	7579312	7579312	C	G	patient_4
chr17	7577127	7577127	C	G	patient_5
chr17	7578190	7578190	T	C	patient_6
chr17	7578534	7578534	C	A	patient_7
chr17	7577106	7577106	G	T	patient_8
chr17	7577022	7577022	G	A	patient_9
```
#### 2. Regions of Interest
Regions of interest can be coding or noncoding should be in a BED12 format

#### 3. Sites of Interest 
Sites of interest (e.g. Transcription factor binding sites) which may be active sites in regions of interest should be in a BED4 format

More thorough documentation of the activePathways function can be found in R with `?activePathways`, and complete tutorials can be found with `browseVignettes(package='ActiveDriverWGSR')`
