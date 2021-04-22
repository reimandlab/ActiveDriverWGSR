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
#'     \item{ref}{the reference allele as a string containing the bases A, T, C, G or -}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C, G or -}
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
#' @param sites A data frame containing the following columns: chr, start, end, id
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the site in base 0 coordinates (BED format)}
#'     \item{end}{the end position of the site in base 0 coordinates (BED format)}
#'     \item{id}{the identifier of the element. id's need to match with those listed in the object elements. }
#' }
#' @param window_size An integer indicating the size of the background window in base pairs that is used to establish
#' the expected mutation rate and respective null model. The default is 50000bps
#' @param filter_hyper_MB Hyper-mutated samples carry many passenger mutations and dilute the signal of true drivers.
#' Samples with a rate greater than \code{filter_hyper_MB} mutations per megabase are excluded.
#' The default is 30 mutations per megabase.
#' @param recovery.dir The directory for storing recovery files. If the directory does not exist, ActiveDriverWGS will create the directory.
#' If the parameter is unspecified, recovery files will not be saved. As an ActiveDriverWGS query for large datasets may be computationally heavy,
#' specifying a recovery directory will recover previously computed results if a query is interrupted.
#'
#' @param mc.cores The number of cores which can be used if multiple cores are available. The default is 1.
#' 
#' @param ref_genome The reference genome used on the analysis. The default option is "hg19", other options are "hg38", "mm9" and "mm10". 
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
#' result = ActiveDriverWGS(mutations = cll_mutations,
#' 		elements = cancer_genes[cancer_genes$id %in% some_genes,])
#' }
ActiveDriverWGS = function(mutations,
							elements,
							sites = NULL,
							window_size = 50000,
							filter_hyper_MB = 30,
							recovery.dir = NULL,
							mc.cores = 1,
							ref_genome = "hg19"){

	# Verifying Format for window_size
	if (!(length(window_size) == 1 && is.numeric(window_size) && window_size > 0)) {
		stop("window size must be a positive integer")
	}
	
	# Verifying Format for mc.cores
	if (!(length(mc.cores) == 1 && is.numeric(mc.cores) && mc.cores > 0)) {
		stop("mc.cores must be a positive integer")
	}
	
	# Verifying Format for filter_hyper_MB
	if (!(length(filter_hyper_MB) == 1 && is.numeric(filter_hyper_MB) && filter_hyper_MB > 0)) {
		stop("filter_hyper_MB must be a positive integer")
	}
	
	# Verifying Format for ref_genome
	permitted_ref_genomes = c("hg19", "hg38", "mm9", "mm10")
	if (is.na(ref_genome) | is.null(ref_genome) | !ref_genome %in% permitted_ref_genomes) {
		stop("ref_genome must one of 'hg19', 'hg38', 'mm9', 'mm10'")
	}
	
	this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
	permitted_chrs = paste0("chr", c(1:22, "X", "Y"))
	
	if (ref_genome == 'hg38') {
		this_genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
		permitted_chrs = paste0("chr", c(1:22, "X", "Y"))
	}
	if (ref_genome == 'mm9') {
		this_genome = BSgenome.Mmusculus.UCSC.mm9::Mmusculus
		permitted_chrs = paste0("chr", c(1:19, "X", "Y"))
	}
	if (ref_genome == 'mm10') {
		this_genome = BSgenome.Mmusculus.UCSC.mm10::Mmusculus
		permitted_chrs = paste0("chr", c(1:19, "X", "Y"))
	}
	
	# Verifying Format for recovery.dir
	if(!is.null(recovery.dir)){
	if(!is.character(recovery.dir) | length(recovery.dir) != 1) stop("recovery.dir must be a string")
	
	if (!dir.exists(recovery.dir)) {
		dir.create(recovery.dir)
			message(paste0("Creating ", recovery.dir))
		}
		
		if(!endsWith(recovery.dir, "[/]") && recovery.dir != ""){
			recovery.dir = paste0(recovery.dir, "/")
		}
	}	
	
	# Verifying Format for Mutations
	if (!is.data.frame(mutations)) {
		stop("mutations must be a data frame")
	}
	
	if (!all(c("chr", "pos1", "pos2", "ref", "alt", "patient") %in% colnames(mutations))) {
		stop("mutations must contain the following columns: chr, pos1, pos2, ref, alt & patient")
	}
	
	if (any(is.na(mutations))) {
		stop("mutations may not contain missing values")
	}
	
	if (any(duplicated(mutations))) {
		stop("duplicated mutations are present. please review your format")
	}
	
	if (!(is.character(mutations$chr) && 
		is.character(mutations$ref) && 
		is.character(mutations$alt))) {
			stop("chr, ref and alt must be character")
	}
	
	if (!all(mutations$chr %in% permitted_chrs)) {
		stop(paste("Only autosomal and sex chromosomes may be used in mutation data (24 for human, 21 for mouse).", 
				"Note that chr23 and chr24 should be formatted as chrX and chrY, respectively"))
	}
	
	if (!(is.numeric(mutations$pos1) && is.numeric(mutations$pos2))) {
	stop("pos1 and pos2 must be numeric")
	}
	
	if (!(all(grepl("[ATGC\\-]", c(mutations$ref, mutations$alt))))) {
	   stop("Reference and alternate alleles must be A, T, C, G or -")
	}
	
	if (!(is.character(mutations$patient))) {
		stop("patient identifier must be a string")
	}
	
	# Creating gr_muts
	mutations = format_muts(mutations = mutations, this_genome = this_genome, filter_hyper_MB = filter_hyper_MB)
	gr_muts = GenomicRanges::GRanges(mutations$chr,
			IRanges::IRanges(mutations$pos1, mutations$pos2), 
			mcols = mutations[,c("patient", "tag")])
	# save(gr_muts, file=paste0(recovery.dir,"gr_muts.rsav"))
	
	# Verifying Format for Elements
	if (!is.data.frame(elements)) {
		stop("elements must be a data frame")
	}
	if (!all(c("chr", "start", "end", "id") %in% colnames(elements))) {
		stop("elements must contain the following columns: chr, start, end & id")
	}
	if (any(is.na(elements))) {
		stop("elements may not contain missing values")
	}
	if (any(duplicated(elements))) {
		stop("duplicated elements are present. please review your format")
	}
	if (!all(elements$chr %in% permitted_chrs)) {
		stop(paste("Only autosomal and sex chromosomes may be used in element coordinates (24 for human, 21 for mouse).", 
				"Note that chr23 and chr24 should be formatted as chrX and chrY, respectively"))
	}
	
	if (!(is.numeric(elements$start) && is.numeric(elements$end))) {
		stop("start and end must be numeric")
	}
	if (!(is.character(elements$id))) {
		stop("element identifier must be a string")
	}
	
	# Creating elements_gr
	gr_element_coords = GenomicRanges::GRanges(elements$chr,
			IRanges::IRanges(elements$start, elements$end),
			mcols = elements[,c("id")])
	
	# Verifying Format for Sites 
	if(!is.null(sites)){
		if (!is.data.frame(sites)) { 
			stop("sites must be a data frame")
		}
		if (!all(c("chr", "start", "end", "id") %in% colnames(sites))) {
			stop("sites must contain the following columns: chr, start, end & id")
		}
		if (any(is.na(sites))) {
			stop("sites may not contain missing values")
		}
		if (any(duplicated(sites))) {
			stop("duplicated sites are present. please review your format")
		}
		if (!all(sites$chr %in% permitted_chrs)) {
			stop(paste("Only autosomal and sex chromosomes may be used in site coordinates (24 for human, 21 for mouse).", 
					"Note that chr23 and chr24 should be formatted as chrX and chrY, respectively"))
		}
		if (!(is.numeric(sites$start) && is.numeric(sites$end))) {
			stop("start and end must be numeric")
		}
		if (!(is.character(sites$id))) {
			stop("site identifier must be a string")
		}
		
		gr_site_coords = GenomicRanges::GRanges(sites$chr,
				IRanges::IRanges(sites$start, sites$end),
				mcols = sites$id)
	} else {
		gr_site_coords = GenomicRanges::GRanges()
	}
	
	# Running ADWGS Test
	all_results = NULL
	
	# Pre-Filtering Results - only elements with 1+ mutation are analyzed, the rest assigned NA
	# for speeding up computation
	muts_elements_overlap = suppressWarnings(GenomicRanges::findOverlaps(gr_element_coords, gr_muts))
	mutated_elements = sort(unique(gr_element_coords$mcols[S4Vectors::queryHits(muts_elements_overlap)]))
	unmutated_elements = sort(unique(gr_element_coords$mcols[!gr_element_coords$mcols %in% mutated_elements]))
	not_done = mutated_elements
	
	# Unmutated Results
	unmutated_results = NULL
	if(length(unmutated_elements) > 0){
		unmutated_results = data.frame(id = unmutated_elements,
				pp_element = NA, element_muts_obs = NA, element_muts_exp = NA, element_enriched = NA,
				pp_site = NA, site_muts_obs = NA, site_muts_exp = NA, site_enriched = NA,
				stringsAsFactors = FALSE)
	}
	if(!(is.null(recovery.dir))){
		unmutated_results$result_number = 1:length(unmutated_elements) + length(mutated_elements)
	}
	cat("Number of Elements with 0 Mutations: ", length(unmutated_elements), "\n")
	
	# Recovered Results
	recovered_results = NULL
	recovered_result_numbers = c()
	if (!is.null(recovery.dir)) {
		results_filenames = list.files(recovery.dir, pattern = "ADWGS_result[0123456789]+_recovery_file.rsav")
		if (length(results_filenames) > 0) {
			results_filenames = paste0("/", results_filenames)
		}
		
		recovered_results = do.call(rbind, lapply(results_filenames, function(filename) {
			load_result = suppressWarnings(try(load(paste0(recovery.dir, filename)), silent = TRUE))
			if (class(load_result) == "try-error") return(NULL)
			result = result
		}))
		recovered_result_numbers = recovered_results$result_number
	}
	
	cat("Tests to do: ", length(not_done), "\n")
	if (length(recovered_result_numbers) > 0) {
		cat("Tests recovered: ", length(unique(recovered_result_numbers)), "\n")
	}
	
	# Mutated Results
	mutated_results = do.call(rbind, parallel::mclapply(1:length(not_done), function(i) {
		if (i %% 100 == 0) {
			cat(i, " elements completed\n")
		}
		
		# skip calculation if this item is completed already
		if (i %in% recovered_result_numbers) {
			return(NULL)
		}
		result = ADWGS_test(
				id = not_done[i], gr_element_coords = gr_element_coords,
				gr_site_coords = gr_site_coords, gr_maf = gr_muts,
				win_size = window_size, this_genome = this_genome)

		# save each result into recovery dir if requested
		if (!is.null(recovery.dir)) {
			result$result_number = i
			save(result, file = paste0(recovery.dir, "ADWGS_result", i, "_recovery_file.rsav"))
		}
		result
	}, mc.cores = mc.cores))
	
	all_results = rbind(recovered_results, mutated_results, unmutated_results)
	
	if (!all.equal(sort(unique(all_results$id)), sort(unique(elements$id)))) {
		stop("Error: Some elements were not evaluated. Results or recovery.dir may be corrupted.\n")
	}
	
	rm(mutated_results, recovered_results, unmutated_results)
	rm(elements, gr_element_coords)
	
	# Formatting Results
	all_results = .fix_all_results(all_results)
	all_results = .get_signf_results(all_results)
	all_results = all_results[,!colnames(all_results) %in% "result_number"]
	all_results
}
