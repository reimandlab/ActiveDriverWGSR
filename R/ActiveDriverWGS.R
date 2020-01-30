# @import GenomicRanges
# @import IRanges
# @import BSgenome.Hsapiens.UCSC.hg19
# @import parallel


#' ActiveDriverWGS is a driver discovery tool for simple somatic mutations in cancer whole genomes
#'
#' @param mutations A data frame containing the following columns: chr, pos1, pos2, ref, alt, patient.
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{pos1}{the start position of the mutation in base 1 coordinates}
#'     \item{pos2}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#'     \item{patient}{the patient identifier as a string}
#' }
#' @param elements A data frame containing the following columns: chr, start, end, id
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the element in base 0 coordinates (BED format)}
#'     \item{end}{the end position of the element in base 0 coordinates (BED format)}
#'     \item{id}{the element identifier - if the element contains multiple segments such as exons,
#'     each segment should be a separate row with the segment coordinates
#'     and the element identifier as id. Elements can be coding or noncoding such as exons of protein
#'     coding genes or active enhancers.}
#' }
#'
#' @param reference The background reference database in interest
#'
#' @param sites A data frame containing the following columns: chr, start, end, id
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the site in base 0 coordinates (BED format)}
#'     \item{end}{the end position of the site in base 0 coordinates (BED format)}
#'     \item{id}{the site identifier - each site should contain only 1 segment and a unique id. If ids are duplicated,
#'     each segment of the site will be treated as an individual site. Sites can be coding or noncoding such as
#'     phosphosites of protein coding genes in genomic coordinates or transcription factor binding sites of active
#'     enhancers.}
#' }
#' @param window_size An integer indicating the size of the background window in base pairs for which the
#' mutation rate is expected to remain the same. The default is 50000bps.
#' @param filter_hyper_MB Hyper-mutated samples carry many passenger mutations and dilute the signal of true drivers.
#' Samples with a rate greater than \code{filter_hyper_MB} mutations per megabase are excluded.
#' The default is 30 mutations per megabase.
#' @param recovery.dir The directory for storing recovery files. If the directory does not exist, ActiveDriverWGS will create the directory.
#' If the parameter is unspecified, recovery files will not be saved. As an ActiveDriverWGS query for large datasets may be computationally heavy,
#' specifying a recovery directory will recover previously computed results if a query is interrupted.
#'
#' @param mc.cores The number of cores which can be used if multiple cores are available. The default is 1.
#'
#' @return A data frame containing the results of driver discovery containing the following columns: id, pp_element,
#' element_muts_obs, element_muts_exp, element_enriched, pp_site, site_muts_obs, site_muts_exp, site_enriched,
#' fdr_element, fdr_site
#' \describe{
#'     \item{id}{A string identifying the element of interest}
#'     \item{pp_element}{The p-value of the element}
#'     \item{element_muts_obs}{The number of patients with a mutation in the element}
#'     \item{element_muts_exp}{The expected number of patients with a mutation in the element with respect to background}
#'     \item{element_enriched}{A boolean indicating whether the element is enriched in mutations}
#'     \item{pp_site}{The p-value of the site}
#'     \item{site_muts_obs}{The number of patients with a mutation in the site}
#'     \item{site_muts_exp}{The expected number of patients with a mutation in the site with respect to element}
#'     \item{site_enriched}{A boolean indicating whether the site is enriched in mutations}
#'     \item{fdr_element}{The FDR corrected p-value of the element}
#'     \item{fdr_site}{The FDR corrected p-value of the site}
#'     \item{has_site_mutations}{A V indicates the presence of site mutations}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(cancer_genes)
#' data(cll_mutations)
#'
#' some_genes = c("ATM", "MYD88", "NOTCH1", "SF3B1", "XPO1",
#' "SOCS1", "CNOT3", "DDX3X", "KMT2A", "HIF1A", "APC")
#'
#' Test with hg19
#'
#' result = ActiveDriverWGS(mutations = cll_mutations,
#' elements = cancer_genes[cancer_genes$id %in% some_genes,], reference = "hg19")
#'
#' Test with hg38
#'
#' some_genestest = paste(some_genes, collapse = ";| ")
#' some_genestest = paste0(some_genestest, ";")
#' elementshg38test = elementshg38[grepl(some_genestest, elementshg38$id),]
#' result = ActiveDriverWGS(mutations = cll_mutations,
#' elements = elementshg38test, reference = "hg38")
#'
#'
#' Test with mm9
#' mouse_genes = c("BRCA1", "BRCA2")
#'
#' mm9cancerelements = elementsmm10[grepl("BRCA1", elementsmm9$id),]
#'
#'
#' Test with mm10
#' mm10cancerelements = elementsmm10[grepl("BRCA1", elementsmm10$id),]
#' }
ActiveDriverWGS = function(mutations,
                           elements,
                           reference,
                           sites = NULL,
                           window_size = 50000,
                           filter_hyper_MB = 30,
                           recovery.dir = NULL,
                           mc.cores = 1){

  # Verifying Format for window_size
  if (!(length(window_size) == 1 && is.numeric(window_size) && window_size > 0)) stop("window size must be a positive integer")

  # Verifying Format for mc.cores
  if (!(length(mc.cores) == 1 && is.numeric(mc.cores) && mc.cores > 0)) stop("mc.cores must be a positive integer")

  # Verifying Format for filter_hyper_MB
  if (!(length(filter_hyper_MB) == 1 && is.numeric(filter_hyper_MB) && filter_hyper_MB > 0)) stop("filter_hyper_MB must be a positive integer")

  # Verifying Format for recovery.dir
  if(!is.null(recovery.dir)){
    if(!is.character(recovery.dir) | length(recovery.dir) != 1) stop("recovery.dir must be a string")

    if (!dir.exists(recovery.dir)){
      dir.create(recovery.dir)
      message(paste0("Creating ", recovery.dir))
    }

    if(!endsWith(recovery.dir, "[/]") && recovery.dir != ""){
      recovery.dir = paste0(recovery.dir, "/")
    }
  } # else{
    # This can be changed to tempdir() if we should need to
    # recovery.dir = "ActiveDriverWGS_recovery"
    # if (!dir.exists(recovery.dir)){
    #   dir.create(recovery.dir)
    # }
    # message(paste0("Writing results to ", recovery.dir))
  # }


  # Checking if reference genome is a valid
  if (!(reference %in% c("hg19", "hg38", "mm9", "mm10"))){stop("reference must be a valid background reference genome")}


  # if reference genome = Homo sapiens, 22 autosomal chromosomes and 2 sex chromosomes
  if (reference == "hg19"){
    chr1 = 1;
    chrY = 24;
    chromosomes <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  }
  if (reference == "hg38"){
    chr1 = 1;
    chrY = 24;
    chromosomes <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  }

  # if reference genome = Mus musculus, 22 autosomal chromosomes and 2 sex chromosomes


  if (reference == "mm9"){
    chr1 = 1;
    chrY = 21;
    chromosomes <- BSgenome.Mmusculus.UCSC.mm9::Mmusculus
  }


  if (reference == "mm10"){
    chr1 = 1;
    chrY = 20;
    chromosomes <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
  }

  # Verifying Format for Mutations
  legal_dna = c('A', 'T', 'C', 'G')
  if (!is.data.frame(mutations)) stop("mutations must be a data frame")
  if (!all(c("chr", "pos1", "pos2", "ref", "alt", "patient") %in% colnames(mutations))) stop("mutations must contain the following columns: chr, pos1, pos2, ref, alt & patient")
  if (any(is.na(mutations))) stop("mutations may not contain missing values")
  if (any(duplicated(mutations))) stop("duplicated mutations are present. please review your format")
  if (!(is.character(mutations$chr) && is.character(mutations$ref) && is.character(mutations$alt))) stop("chr, ref and alt must be character")
  if (!any(mutations$chr %in% BSgenome::seqnames(chromosomes)[chr1:chrY]))
    stop("Only the 22 autosomal and 2 sex chromosomes may be used at this time.
         Note that chr23 should be formatted as chrX and chr24 should be formatted as chrY")
  if (!(is.numeric(mutations$pos1) && is.numeric(mutations$pos2))) stop("pos1 and pos2 must be numeric")
  if (!(grepl("^[ATCG]+$", mutations$ref) && grepl("^[ATCG]+$", mutations$alt))) stop("Reference and alternate alleles must be A, T, C or G")
  if (!(is.character(mutations$patient))) stop("patient identifier must be a string")

  # Creating gr_muts
  mutations = format_muts(mutations = mutations,
                          filter_hyper_MB = filter_hyper_MB,
                          reference = chromosomes)
  gr_muts = GenomicRanges::GRanges(mutations$chr,
                                   IRanges::IRanges(mutations$pos1, mutations$pos2), mcols=mutations[,c("patient", "tag")])
  # save(gr_muts, file=paste0(recovery.dir,"gr_muts.rsav"))

  # Verifying Format for Elements
  if (!is.data.frame(elements)) stop("elements must be a data frame")
  if (!all(c("chr", "start", "end", "id") %in% colnames(elements))) stop("elements must contain the following columns: chr, start, end & id")
  if (any(is.na(elements))) stop("elements may not contain missing values")
  if (any(duplicated(elements))) stop("duplicated elements are present. please review your format")
  # if (!all(elements$chr %in% BSgenome.Hsapiens.UCSC.hg19::seqnames(Hsapiens)[1:24])) stop("Only the 22 autosomal and 2 sex chromosomes may be used at this time. Note that chr23 should be formatted as chrX and chr24 should be formatted as chrY")
  if (!(is.numeric(elements$start) && is.numeric(elements$end))) stop("start and end must be numeric")
  if (!(is.character(elements$id))) stop("element identifier must be a string")

  # Creating elements_gr
  gr_element_coords = GenomicRanges::GRanges(elements$chr,
                                             IRanges::IRanges(elements$start, elements$end),
                                             mcols = elements[,c("id")])

  # Verifying Format for Sites
  if(!is.null(sites)){
    if (!is.data.frame(sites)) stop("sites must be a data frame")
    if (!all(c("chr", "start", "end", "id") %in% colnames(sites))) stop("sites must contain the following columns: chr, start, end & id")
    if (any(is.na(sites))) stop("sites may not contain missing values")
    if (any(duplicated(sites))) stop("duplicated sites are present. please review your format")
    # if (!all(sites$chr %in% BSgenome.Hsapiens.UCSC.hg19::seqnames(Hsapiens)[1:24])) stop("Only the 22 autosomal and 2 sex chromosomes may be used at this time. Note that chr23 should be formatted as chrX and chr24 should be formatted as chrY")
    if (!(is.numeric(sites$start) && is.numeric(sites$end))) stop("start and end must be numeric")
    if (!(is.character(sites$id))) stop("site identifier must be a string")

    gr_site_coords = GenomicRanges::GRanges(sites$chr,
                                            IRanges::IRanges(sites$start, sites$end))
  }else{
    gr_site_coords = GenomicRanges::GRanges()
    # gr_site_coords = GenomicRanges::GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
  }

  # Running ADWGS Test
  all_results = NULL
  not_done = sort(unique(gr_element_coords$mcols))

  recovered_results = NULL
  recovered_result_numbers = c()
  if (!is.null(recovery.dir)) {
    results_filenames = list.files(recovery.dir, pattern = "ADWGS_result[0123456789]+_recovery_file.rsav")
    if (length(results_filenames) > 0) results_filenames = paste0("/", results_filenames)
    recovered_results = do.call(rbind, lapply(results_filenames, function(filename) {
      load_result = suppressWarnings(try(load(paste0(recovery.dir, filename)), silent = TRUE))
      if (class(load_result) == "try-error") return(NULL)
      if (ncol(result) != 13) return(NULL)
      result = result
    }))
    recovered_result_numbers = recovered_results$result_number
  }

  cat("Tests to do: ", length(not_done), "\n")
  if (length(recovered_result_numbers) > 0) cat("Tests recovered: ", length(recovered_result_numbers), "\n")

  mutated_results = do.call(rbind, parallel::mclapply(1:length(not_done), function(i) {
    if (i %% 100 == 0) cat(i, " elements completed\n")
    if (i %in% recovered_result_numbers) return(NULL)
    result = ADWGS_test(id = not_done[i],
                        gr_element_coords = gr_element_coords,
                        gr_site_coords = gr_site_coords,
                        gr_maf = gr_muts,
                        win_size = window_size,
                        background = chromosomes)
    if (!is.null(recovery.dir)) {
      result$result_number = i
      save(result, file = paste0(recovery.dir, "ADWGS_result", i, "_recovery_file.rsav"))
    }
    result
  }, mc.cores = mc.cores))

  all_results = rbind(recovered_results, mutated_results)

  rm(mutated_results, recovered_results)
  if (nrow(all_results) != length(unique(elements$id))) stop("Error: Something unexpected happened. Please try again.\n")
  rm(elements, gr_element_coords)

  # Formatting Results
  all_results = .fix_all_results(all_results)
  all_results = .get_signf_results(all_results)
  all_results = all_results[,!colnames(all_results) %in% "result_number"]
  all_results
}
