## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      warning=FALSE, 
                      message=FALSE, 
                      width=500)
options(max.print=35)

## ----pressure, echo=FALSE, fig.cap="The ActiveDriverWGS Model", out.width = '75%'----
knitr::include_graphics("ADWGS_diagram.png")

## ----input--------------------------------------------------------------------
library(ActiveDriverWGS)

data("cll_mutations")
head(cll_mutations)

data("cancer_genes")
head(cancer_genes)

data("cancer_gene_sites")
head(cancer_gene_sites)


## ----prepare_elements_from_BED12 and prepare_elements_from_BED4---------------
elements = prepare_elements_from_BED12(
  system.file(
    "extdata", 
    "chr17.coding_regions.bed", 
    package = "ActiveDriverWGS", 
    mustWork = TRUE))

head(elements)

sites = prepare_elements_from_BED4(
  system.file(
    "extdata", 
    "chr17.PTM_sites.bed", 
    package = "ActiveDriverWGS", 
    mustWork = TRUE))

head(sites)

## ----ActiveDriverWGS----------------------------------------------------------

some_genes = c("ATM", "MYD88", "NOTCH1")

results = ActiveDriverWGS(mutations = cll_mutations,
                          elements = cancer_genes[cancer_genes$id %in% some_genes,],
                          sites = cancer_gene_sites)


## ----results------------------------------------------------------------------
results


## ----pipeline-----------------------------------------------------------------
library(GenomicRanges)

# Loading elements & creating a GRanges object
data(cancer_genes)
gr_element_coords = GRanges(seqnames = cancer_genes$chr,
                            IRanges(start = cancer_genes$start,
                                    end = cancer_genes$end),
                            mcols = cancer_genes$id)

# Loading sites & creating a GRanges object
data(cancer_gene_sites)
gr_site_coords = GRanges(seqnames = cancer_gene_sites$chr,
                         IRanges(start = cancer_gene_sites$start,
                                 end = cancer_gene_sites$end),
                         mocols = cancer_gene_sites$id)

# Loading mutations, format muts & creating a GRanges object
data(cll_mutations)

# load the default reference genome
this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens

# format_muts
cll_mutations = format_muts(cll_mutations, this_genome,
                            filter_hyper_MB = 30)

gr_maf = GRanges(cll_mutations$chr,
                 IRanges(start = cll_mutations$pos1,
                         end = cll_mutations$pos2),
                 mcols = cll_mutations[,c("patient", "tag")])

# Examplifying the ATM Element
id = "ATM"

## ----result-------------------------------------------------------------------
# Result of 1 input element
result = ADWGS_test(id = id,
                    gr_element_coords = gr_element_coords,
                    gr_site_coords = gr_site_coords,
                    gr_maf = gr_maf,
                    win_size = 50000, this_genome = this_genome)

result

