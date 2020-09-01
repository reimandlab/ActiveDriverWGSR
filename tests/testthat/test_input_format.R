context("Format of Input Files to ActiveDriverWGS")

set.seed(100)

# loading regions
data(cancer_genes)

# loading sites
data(cancer_gene_sites)

# loading mutations
data(cll_mutations)
cll_mutations = cll_mutations[sample(nrow(cll_mutations), 30),]

# Some CLL drivers in this dataset + random genes
some_genes = c("ATM", "MYD88",
               "NOTCH1","SF3B1",
               "XPO1", "SOCS1",
               "CNOT3", "DDX3X",
               "KMT2A", "HIF1A",
               "APC")

test_that("mutations are formatted properly",{

  # Mutations are in a dataframe
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = as.matrix(cll_mutations)),
               "mutations must be a data frame")

  # Colnames of the mutations
  this_mutations = cll_mutations
  colnames(this_mutations) = paste0("V", 1:6)
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "mutations must contain the following columns: chr, pos1, pos2, ref, alt & patient")

  # Colnames of the mutations
  this_mutations = cll_mutations[,c(1:2,4:6)]
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "mutations must contain the following columns: chr, pos1, pos2, ref, alt & patient")

  # Missing values
  this_mutations = cll_mutations
  this_mutations[sample(1:nrow(this_mutations), 10),3] = NA
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "mutations may not contain missing values")

  # Duplicated mutations
  this_mutations = cll_mutations
  this_mutations = rbind(this_mutations, this_mutations[sample(1:nrow(this_mutations), 10),])
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "duplicated mutations are present. please review your format")

  # Illegal chromosomes
  this_chromosomes = structure(paste0("chr", letters[1:24]), names = paste0("chr", c(1:22, "X", "Y")))
  this_mutations = cll_mutations
  this_mutations$chr = this_chromosomes[this_mutations$chr]
  err_expect = paste("Only autosomal and sex chromosomes may be used in mutation data (24 for human, 21 for mouse).", 
				"Note that chr23 and chr24 should be formatted as chrX and chrY, respectively")
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               err_expect, fixed = TRUE)

  # Illegal bases
  this_mutations = cll_mutations
  this_mutations$ref = tolower(this_mutations$ref)
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "Reference and alternate alleles must be A, T, C, G or -")

  this_mutations = cll_mutations
  this_mutations$alt = tolower(this_mutations$alt)
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "Reference and alternate alleles must be A, T, C, G or -")

  # Positions as strings
  this_mutations = cll_mutations
  this_mutations$pos1 = as.character(this_mutations$pos1)
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "pos1 and pos2 must be numeric")

  this_mutations = cll_mutations
  this_mutations$pos2 = as.character(this_mutations$pos2)
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "pos1 and pos2 must be numeric")

  # Patient identifier as numeric
  this_num_patient = structure(1:length(unique(cll_mutations$patient)), names = unique(cll_mutations$patient))
  this_mutations = cll_mutations
  this_mutations$patient = as.numeric(this_num_patient[this_mutations$patient])
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "patient identifier must be a string")

  # Chromosomes are numeric
  this_num_chrom = structure(1:24, names = paste0("chr", c(1:22, "X", "Y")))
  this_mutations = cll_mutations
  this_mutations$chr = as.numeric(this_num_chrom[this_mutations$chr])
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = this_mutations),
               "chr, ref and alt must be character")

})

test_that("elements are formatted properly",{

  # Elements are in a dataframe
  expect_error(ActiveDriverWGS(elements = as.matrix(cancer_genes),
                               mutations = cll_mutations),
               "elements must be a data frame")

  # Elements contain the right columns
  expect_error(ActiveDriverWGS(elements = cancer_genes[,1:3],
                               mutations = cll_mutations),
               "elements must contain the following columns: chr, start, end & id")

  # Elements does not contain NAs
  this_elements = cancer_genes
  this_elements[1,1] = NA
  expect_error(ActiveDriverWGS(elements = this_elements,
                               mutations = cll_mutations),
               "elements may not contain missing values")

  # Duplicated elements
  this_elements = rbind(cancer_genes[1:10,], cancer_genes[1:10,])
  expect_error(ActiveDriverWGS(elements = this_elements,
                               mutations = cll_mutations),
               "duplicated elements are present. please review your format")

  # Start End non-numeric
  this_elements = cancer_genes
  this_elements$start = as.character(this_elements$start)
  expect_error(ActiveDriverWGS(elements = this_elements,
                               mutations = cll_mutations),
               "start and end must be numeric")

  # unexpected chromosome
  this_elements = cancer_genes
  this_elements$chr[1] = paste0("CHR", this_elements$chr[1])
  err_expect = paste("Only autosomal and sex chromosomes may be used in element coordinates (24 for human, 21 for mouse).", 
			"Note that chr23 and chr24 should be formatted as chrX and chrY, respectively")
  expect_error(ActiveDriverWGS(elements = this_elements,
                               mutations = cll_mutations),
	err_expect, fixed = TRUE)

  # Numeric ID
  this_elements = cancer_genes
  this_elements$id = 1:nrow(this_elements)
  expect_error(ActiveDriverWGS(elements = this_elements,
                               mutations = cll_mutations),
               "element identifier must be a string")
})

test_that("sites are formatted properly",{

  # Sites are in a dataframe
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               sites = as.matrix(cancer_gene_sites)),
               "sites must be a data frame")

  # Sites contain the right columns
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               sites = cancer_gene_sites[,1:3]),
               "sites must contain the following columns: chr, start, end & id")


  # Sites does not contain NAs
  this_sites = cancer_gene_sites
  this_sites[1,1] = NA
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               sites = this_sites),
               "sites may not contain missing values")

  # Duplicated sites
  this_sites = rbind(cancer_gene_sites[1:10,], cancer_gene_sites[1:10,])
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               sites = this_sites),
               "duplicated sites are present. please review your format")

  # Start End non-numeric
  this_sites = cancer_gene_sites
  this_sites$start = as.character(this_sites$start)
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               sites = this_sites),
               "start and end must be numeric")

  # Numeric ID
  this_sites = cancer_gene_sites
  this_sites$id = 1:nrow(this_sites)
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               sites = this_sites),
               "site identifier must be a string")

  # unexpected chromosome
  this_sites = cancer_gene_sites
  this_sites$chr[1] = paste0("CHR", this_sites$chr[1])
  err_expect = paste("Only autosomal and sex chromosomes may be used in site coordinates (24 for human, 21 for mouse).", 
		"Note that chr23 and chr24 should be formatted as chrX and chrY, respectively")
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               sites = this_sites),
		err_expect, fixed = TRUE)
})

test_that("window_size is a positive integer",{

  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               window_size = -10000),
               "window size must be a positive integer")

  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               window_size = "50000"),
               "window size must be a positive integer")

})

test_that("filter_hyper_MB is a positive integer",{

  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               filter_hyper_MB = -2),
               "filter_hyper_MB must be a positive integer")

  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               filter_hyper_MB = "30"),
               "filter_hyper_MB must be a positive integer")

})

test_that("recover.dir is appropriate",{
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               recovery.dir = 2),
               "recovery.dir must be a string")
})

test_that("mc.cores is a positive integer",{
  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               mc.cores = -2),
               "mc.cores must be a positive integer")

  expect_error(ActiveDriverWGS(elements = cancer_genes,
                               mutations = cll_mutations,
                               mc.cores = "a"),
               "mc.cores must be a positive integer")
})

