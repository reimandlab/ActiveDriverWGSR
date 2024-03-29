context("Testing the function of the format_muts function")

# loading mutations
data(cll_mutations)

# Formatting muts
this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
formatted_muts = format_muts(cll_mutations[1:10,], this_genome = this_genome)

test_that("format_muts returns a data frame with the right columns",{

  this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  expect_warning(format_muts(cll_mutations, this_genome = this_genome), "reversing 0 positions")

  expect_identical(colnames(formatted_muts)[ncol(formatted_muts)], "tag")

  expect_is(formatted_muts[,"chr"], "character")

  expect_is(formatted_muts[,"pos1"], "numeric")

  expect_is(formatted_muts[,"pos2"], "numeric")

  expect_is(formatted_muts[,"ref"], "character")

  expect_is(formatted_muts[,"alt"], "character")

  expect_is(formatted_muts[,"patient"], "character")

  expect_is(formatted_muts[,"tag"], "character")

  # Testing that filtering works
  some_patients = c("001-0002-03TD", "003-0005-09TD", "012-02-1TD", "125", "128", "141", "178")
  this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  expect_warning(format_muts(cll_mutations[cll_mutations$patient %in% some_patients,], 
  							this_genome = this_genome, filter_hyper_MB = 1),
                "remove 2 hypermut samples, n=6709 muts, 50%")

  # Testing that SNVs only works
  this_mutations = cll_mutations[cll_mutations$pos2 == cll_mutations$pos1,]
  this_mutations = this_mutations[1:10,]
  this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  test_snvs_only = format_muts(mutations = this_mutations, this_genome = this_genome)
  expect_identical(colnames(test_snvs_only)[ncol(test_snvs_only)], "tag")
})

test_that("testing errors on the format muts function",{

  # Filtering hypermutated samples
  this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  expect_error(format_muts(mutations = cll_mutations, this_genome = this_genome,
                           filter_hyper_MB = 0.1),
                 "No mutations left after filtering hypermutators")

  # Testing that mutations in unsequenceable regions are filtered (Part of .get_3n_context_of_mutations)
  this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  this_mutations = data.frame("chr" = "chr1",
                              "pos1" = 126000000,
                              "pos2" = 126000000,
                              "ref" = "C",
                              "alt" = "A",
                              "patient" = "Lady Gaga",
                              stringsAsFactors = F)
  this_mutations = rbind(cll_mutations[1:10,], this_mutations)
  expect_warning(format_muts(mutations = this_mutations, this_genome = this_genome,
                            filter_hyper_MB = 30),
                "Removing 1 invalid SNVs & indels")

  # This test doesn't work - I'm not sure why yet
  # Testing that mutations outside of ranges do not work (Part of .get_3n_context_of_mutations)
  # this_mutations = data.frame("chr" = "chr1",
  #                             "pos1" = 249250621,
  #                             "pos2" = 249250622,
  #                             "ref" = "C",
  #                             "alt" = "A",
  #                             "patient" = "Lady Gaga",
  #                             stringsAsFactors = F)
  # expect_error(format_muts(mutations = this_mutations,
  #                          filter_hyper_MB = 30), "")

})
