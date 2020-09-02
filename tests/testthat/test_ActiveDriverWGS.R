context("Testing the ActiveDriverWGS function")

library(GenomicRanges)

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
this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
cll_mutations = format_muts(cll_mutations, this_genome = this_genome,
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

results = ActiveDriverWGS(elements = cancer_genes[cancer_genes$id %in% some_genes,],
                          mutations = cll_mutations,
                          sites = cancer_gene_sites)

test_that("Results have the appropriate format",{
  
  # results are a data.frame
  expect_is(results, "data.frame")
  
  # results have the write column names
  expect_identical(colnames(results), c("id", "pp_element", "element_muts_obs", "element_muts_exp", "element_enriched", "pp_site",           
                                        "site_muts_obs", "site_muts_exp", "site_enriched", "fdr_element", "fdr_site", "has_site_mutations"))
  
  # id is character
  expect_is(results[,"id"], "character")
  
  # pp_element is numeric
  expect_is(results[,"pp_element"], "numeric")
  
  # pp_element  < 1
  expect_true(all(results[,"pp_element"] <= 1))
  
  # element_muts_obs is numeric
  expect_is(results[,"element_muts_obs"], "numeric")
  
  # element_muts_obs has the right number of patients
  for(i in 1:nrow(results)){
    num_patients = length(unique(gr_maf[queryHits(findOverlaps(gr_maf, gr_element_coords[gr_element_coords$cols == results$id[i]]))]$mcols.patient))
    expect_equal(results$element_muts_obs[i], num_patients)
  }

  # element_muts_exp is numeric
  expect_is(results[,"element_muts_exp"], "numeric")
  
  # element_enriched
  expect_is(results[,"element_enriched"], "logical")
  
  # pp_site is numeric
  expect_is(results[,"pp_site"], "numeric")
  
  # pp_site < 1
  expect_true(all(results[,"pp_site"] <= 1))
  
  # site_enriched
  expect_is(results[,"site_enriched"], "logical")
  
  # fdr_element
  expect_is(results[,"fdr_element"], "numeric")
  
  # fdr_element < 1
  expect_true(all(results[,"fdr_element"] <= 1))
  
  # fdr_site
  expect_is(results[,"fdr_site"], "numeric")
  
  # fdr_site < 1
  expect_true(all(results[,"fdr_site"] <= 1))
  
  # has_site_mutations
  expect_is(results[,"has_site_mutations"], "character")
  

  
})
