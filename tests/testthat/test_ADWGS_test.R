context("Testing the ADWGS_test function")

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

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
                         mcols = gsub("(.+) (.+) (.+)", "\\1", cancer_gene_sites$id))
# loading mutations
data(cll_mutations)

this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
cll_mutations = format_muts(cll_mutations, this_genome = this_genome,
                            filter_hyper_MB = 30)

gr_maf = GRanges(cll_mutations$chr,
                 IRanges(start = cll_mutations$pos1,
                         end = cll_mutations$pos2),
                 mcols=cll_mutations[,c("patient", "tag")])

id = "XPO1"

result = ADWGS_test(id, gr_element_coords, gr_site_coords, gr_maf, 50000, this_genome)

test_that("filter_hyper_MB is a positive integer",{

  # result are a data.frame
  expect_is(result, "data.frame")

  # result have the write column names
  expect_identical(colnames(result), c("id", "pp_element", "element_muts_obs", "element_muts_exp", "element_enriched", "pp_site",
                                        "site_muts_obs", "site_muts_exp", "site_enriched"))

  # id is character
  expect_is(result[,"id"], "character")

  # pp_element is numeric
  expect_is(result[,"pp_element"], "numeric")

  # pp_element  < 1
  expect_true(all(result[,"pp_element"] <= 1))

  # element_muts_obs is numeric
  expect_is(result[,"element_muts_obs"], "numeric")

  # element_muts_obs has the right number of patients
  num_patients = length(unique(gr_maf[queryHits(findOverlaps(gr_maf, gr_element_coords[gr_element_coords$cols == result$id[1]]))]$mcols.patient))
  expect_equal(result$element_muts_obs[1], num_patients)

  # element_muts_exp is numeric
  expect_is(result[,"element_muts_exp"], "numeric")

  # element_enriched
  expect_is(result[,"element_enriched"], "logical")

  # pp_site is numeric
  expect_is(result[,"pp_site"], "numeric")

  # pp_site < 1
  expect_true(all(result[,"pp_site"] <= 1))

  # site_enriched
  expect_is(result[,"site_enriched"], "logical")


})
