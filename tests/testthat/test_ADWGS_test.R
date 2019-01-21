context("Testing the ADWGS_test function")

library(ActiveDriverWGSR)
library(testthat)

# loading regions
data(cancer_genes)
gr_element_coords = GRanges(seqnames = cancer_genes$chr, 
                            IRanges(start = cancer_genes$start, 
                                    end = cancer_genes$end),
                            cols = cancer_genes$id)

# loading sites
data(cancer_gene_sites)
gr_site_coords = GRanges(seqnames = cancer_gene_sites$chr,
                         IRanges(start = cancer_gene_sites$start,
                                 end = cancer_gene_sites$end),
                         mocols = cancer_gene_sites$id)
# loading mutations
data(cll_mutations)

cll_mutations = format_muts(cll_mutations,
                            filter_hyper_MB = 30)

gr_maf = GRanges(cll_mutations$chr,
                 IRanges(start = cll_mutations$pos1,
                         end = cll_mutations$pos2),
                 mcols=cll_mutations[,c("patient", "tag")])

# Some CLL drivers in this dataset + random genes
some_genes = c("ATM", "MYD88",
               "NOTCH1","SF3B1",
               "XPO1", "SOCS1",
               "CNOT3", "DDX3X",
               "KMT2A", "HIF1A",
               "APC")

id = "ATM"

test_that("filter_hyper_MB is a positive integer",{
  # testing ADWGS_test
  
})
