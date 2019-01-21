context("Testing the function of the format_muts function")

library(ActiveDriverWGSR)
library(testthat)

# loading mutations
data(cll_mutations)

# Some CLL drivers in this dataset + random genes
some_genes = c("ATM", "MYD88",
               "NOTCH1","SF3B1",
               "XPO1", "SOCS1",
               "CNOT3", "DDX3X",
               "KMT2A", "HIF1A",
               "APC")


test_that("testing the format_muts function",{
  
  # Filtering hypermutated samples
  expect_error(format_muts(mutations = cll_mutations,
                           filter_hyper_MB = 0.1),
                 "No mutations left after filtering")
  
  # Illegal mutations
  this_mutations = cll_mutations
  this_mutations$chr = gsub("chr", "ABC", this_mutations$chr)
  expect_error(format_muts(mutations = this_mutations,
                           filter_hyper_MB = 30),
               "No mutations left after filtering")
  
  # Testing get_3n_context_of_mutations
  
})
