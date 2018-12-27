require(parallel)

#' ActiveDriverWGS is a driver discovery tool for simple somatic mutations in cancer whole genomes
#'
#' @param mutations 
#' @param elements 
#' @param sites 
#' @param window_size 
#' @param filter_hyper_MB 
#' @param recovery.dir 
#' @param mc.cores 
#'
#' @return
#' @export
#'
#' @examples
#' TODO Make sure all mutations are legal
#' TODO Format all results
ActiveDriverWGS = function(mutations, 
                           elements, 
                           sites = NULL, 
                           window_size = 50000,
                           filter_hyper_MB = 30,
                           recovery.dir = "",
                           mc.cores = 1){
  
  # Verifying Format for window_size
  if (!(length(window_size) == 1 && is.numeric(window_size) && window_size > 0)) stop("window size must be a positive integer")

  # Verifying Format for mc.cores
  if (!(length(mc.cores) == 1 && is.numeric(mc.cores) && mc.cores > 0)) stop("mc.cores must be a positive integer")
  
  # Verifying Format for filter_hyper_MB
  if (!(length(filter_hyper_MB) == 1 && is.numeric(filter_hyper_MB) && filter_hyper_MB > 0)) stop("filter_hyper_MB must be a positive integer")
  
  # Verifying Format for recovery.dir
  if(!is.character(recovery.dir) | length(recovery.dir) != 1) stop("recovery.dir must be a string")
  if(recovery.dir != ""){
    if (!dir.exists(recovery.dir)){
      dir.create(recovery.dir)
      message(paste0("Creating ", recovery.dir))
    }
  }
  if(!endsWith(recovery.dir, "[/]")){
    recovery.dir = paste0(recovery.dir, "/")
  }
  
  # Verifying Format for Mutations
  legal_dna = c('A', 'T', 'C', 'G')
  if (!is.data.frame(mutations)) stop("mutations must be a data frame")
  if (!all(c("chr", "pos1", "pos2", "ref", "alt", "patient") %in% colnames(mutations))) stop("mutations must contain the following columns: chr, pos1, pos2, ref, alt & patient")
  if (any(is.na(mutations))) stop("mutations may not contain missing values")
  if (any(duplicated(mutations))) stop("duplicated mutations are present. please review your format")
  if (!all(mutations$chr %in% seqnames(Hsapiens)[1:24])) stop("Only the 22 autosomal and 2 sex chromosomes may be used at this time. Note that chr23 should be formatted as chrX and chr24 should be formatted as chrY")
  if (!(is.integer(mutations$pos1) && is.integer(mutations$pos2))) stop("pos1 and pos2 must be integers")
  if (!(all(mutations$ref %in% legal_dna) && all(mutations$alt %in% legal_dna))) stop("Reference and alternate alleles must be A, T, C or G")
  if (!(is.character(mutations$patient))) stop("patient identifier must be a string")
  
  # Creating gr_muts
  mutations = format_muts(mutations, filter_hyper_MB = filter_hyper_MB)
  gr_muts = GRanges(mutations$chr, 
                    IRanges(mutations$pos1, mutations$pos2), mcols=mutations[,c("patient", "tag")])
  save(gr_muts, file=paste0(recovery.dir,"/gr_muts.rsav"))
  
  # Verifying Format for Elements
  if (!is.data.frame(elements)) stop("elements must be a data frame")
  if (!all(c("chr", "start", "end", "id") %in% colnames(elements))) stop("elements must contain the following columns: chr, start, end & id")
  if (any(is.na(elements))) stop("elements may not contain missing values")
  if (any(duplicated(elements))) stop("duplicated elements are present. please review your format")
  if (!all(elements$chr %in% seqnames(Hsapiens)[1:24])) stop("Only the 22 autosomal and 2 sex chromosomes may be used at this time. Note that chr23 should be formatted as chrX and chr24 should be formatted as chrY")
  if (!(is.integer(elements$start) && is.integer(elements$end))) stop("start and end must be integers")
  if (!(is.character(elements$id))) stop("element identifier must be a string")
  
  # Creating elements_gr
  gr_element_coords = GRanges(elements$chr,
                              IRanges(elements$start, elements$end),
                              mcols = elements[,c("id")])
  
  # Verifying Format for Sites
  if(!is.null(sites)){
    if (!is.data.frame(sites)) stop("sites must be a data frame")
    if (!all(c("chr", "start", "end", "id") %in% colnames(sites))) stop("sites must contain the following columns: chr, start, end & id")
    if (any(is.na(sites))) stop("sites may not contain missing values")
    if (any(duplicated(sites))) stop("duplicated sites are present. please review your format")
    if (!all(sites$chr %in% seqnames(Hsapiens)[1:24])) stop("Only the 22 autosomal and 2 sex chromosomes may be used at this time. Note that chr23 should be formatted as chrX and chr24 should be formatted as chrY")
    if (!(is.integer(sites$start) && is.integer(sites$end))) stop("start and end must be integers")
    if (!(is.character(sites$id))) stop("site identifier must be a string")
    
    gr_site_coords = GRanges(sites$chr,
                             IRanges(sites$start, sites$end))
  }else{
    gr_site_coords = GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
  }
  
  # Running ADWGS Test
  all_results = NULL
  not_done = sort(unique(gr_element_coords$mcols))
  
  recovered_results = NULL
  recovered_result_numbers = c()
  if (!is.na(recovery_dir)) {
    results_filenames = list.files(recovery_dir, pattern = "ADWGS_result[0123456789]+_recovery_file.rsav")
    if (length(results_filenames) > 0) results_filenames = paste0("/", results_filenames)
    recovered_results = do.call(rbind, lapply(results_filenames, function(filename) {
      load_result = suppressWarnings(try(load(paste0(recovery_dir, filename)), silent = TRUE))
      if (class(load_result) == try_error) return(NULL)
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
    element_id = not_done[i]
    result = ADWGS_test(id = not_done[i],
                        gr_element_coords = gr_element_coords,
                        gr_site_coords = gr_site_coords,
                        gr_maf = gr_muts,
                        win_size = window_size)
    if (!is.na(recovery_dir)) {
      result$result_number = i
      save(result, file = paste0(recovery_dir, "/ADWGS_result", i, "_recovery_file.rsav"))
    }
    result
  }, mc.cores = mc.cores))
  
  mutated_results = rbind(recovered_results, mutated_results)
  
  # Calculating FDR Element, FDR Site and Ordering from Most to Least Significant
  
  mutated_results
}