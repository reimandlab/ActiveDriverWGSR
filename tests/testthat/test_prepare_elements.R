context("Testing the prepare_elements_from_BED12 function")

test_that("Prepare Elements Imports Elements in the Correct Format", {

  elements = prepare_elements_from_BED12(system.file("extdata",
                                                     "chr17.coding_regions.bed",
                                                     package = "ActiveDriverWGSR",
                                                     mustWork = TRUE))

  expect_is(elements, "data.frame")

  expect_equal(ncol(elements), 4)

  expect_is(elements[,1], "character")

  expect_is(elements[,2], "numeric")

  expect_is(elements[,3], "numeric")

  expect_is(elements[,4], "character")

  expect_identical(colnames(elements), c("chr", "starts", "ends", "id"))

})

test_that("Prepare Elements Throws the Right Errors",{

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
