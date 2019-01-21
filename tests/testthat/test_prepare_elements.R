context("Testing the prepare_elements_from_BED12 function")

library(ActiveDriverWGSR)
library(testthat)

test_that("filter_hyper_MB is a positive integer",{
  
  expect_error(prepare_elements_from_BED12(
    system.file(
      "extdata", 
      "cds_column_error.bed", 
      package = "ActiveDriverWGSR", 
      mustWork = TRUE)),
      "Incorrect BED12 Format: 12 Columns in BED12 files")
  
  expect_error(prepare_elements_from_BED12(
    system.file(
      "extdata", 
      "cds_str_coords_error.bed", 
      package = "ActiveDriverWGSR", 
      mustWork = TRUE)),
    "Incorrect BED12 Format: Incorrect coordinate format")
  
  expect_error(prepare_elements_from_BED12(
    system.file(
      "extdata", 
      "cds_chrom_error.bed", 
      package = "ActiveDriverWGSR", 
      mustWork = TRUE)),
    "Incorrect BED12 Format: Chromosomes must be autosomal, sex or mitochondrial")
  
  expect_error(prepare_elements_from_BED12(
    system.file(
      "extdata", 
      "cds_id_error.bed", 
      package = "ActiveDriverWGSR", 
      mustWork = TRUE)),
    "Incorrect BED12 Format: IDs must be a character string")
  
  expect_error(prepare_elements_from_BED12(
    system.file(
      "extdata", 
      "cds_blockcount_error.bed", 
      package = "ActiveDriverWGSR", 
      mustWork = TRUE)),
    "Incorrect BED12 Format: Incorrect blockCounts")
  
  expect_error(prepare_elements_from_BED12(
    system.file(
      "extdata", 
      "cds_blockstart_error.bed", 
      package = "ActiveDriverWGSR", 
      mustWork = TRUE)),
    "Incorrect BED12 Format: Incorrect blockSizes or blockStarts")
  
  expect_error(prepare_elements_from_BED12(
    system.file(
      "extdata", 
      "cds_blocksize_error.bed", 
      package = "ActiveDriverWGSR", 
      mustWork = TRUE)),
    "Incorrect BED12 Format: Incorrect blockSizes or blockStarts")

})
