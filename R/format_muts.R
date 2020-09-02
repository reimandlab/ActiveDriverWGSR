

# @import GenomicRanges
# @import IRanges
# @import BSgenome
# @import Biostrings

#' This function finds the tri-nucleotide context of mutations
#'
#' @param mutations A data frame with the following columns: chr, pos1, pos2, ref, alt, patient
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{pos1}{the start position of the mutation in base 1 coordinates}
#'     \item{pos2}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#'     \item{patient}{the patient identifier as a string}
#' }
#' @param this_genome The reference genome object of BSgenome, for example BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#'
#' @return A data frame consisting of the same columns as the original mutations data frame and sorted
#' by SNVs and Indels with an additional column \code{tag} which indicates the trinucleotide context of
#' the mutation
.get_3n_context_of_mutations = function(mutations, this_genome) {
	
	legal_dna = c("A", "C", "G", "T")
	
	mutations_snv = mutations[mutations$ref %in% legal_dna & mutations$alt %in% legal_dna,]
	mutations_indel = mutations[!(mutations$ref %in% legal_dna & mutations$alt %in% legal_dna),]
	
	if (nrow(mutations_snv) > 0) {
		# snvs can have flanks
		flank_ranges = GenomicRanges::GRanges(mutations_snv$chr,
				IRanges::IRanges(start = mutations_snv$pos1-1, 
								end = mutations_snv$pos2+1), strand="*")
		triples = as.character(BSgenome::getSeq(this_genome, flank_ranges))
		
		# filtering SNVs in unsequenceable regions of the genome
		seq_snvs = grep("N", triples, invert = T)
		mutations_snv = mutations_snv[seq_snvs,]
		triples = triples[seq_snvs]
		
		# complement trinucleotide where necessary to force into one signature space
		new_triples = triples
		new_alt = mutations_snv$alt
		which_to_complement = which(mutations_snv$ref %in% c("G", "A"))
		
		new_triples[which_to_complement] = 
				as.character(Biostrings::complement(Biostrings::DNAStringSet(
					new_triples[which_to_complement])))
		new_alt[which_to_complement] = 
				as.character(Biostrings::complement(Biostrings::DNAStringSet(
					new_alt[which_to_complement])))
		mutations_snv$tag = paste0(new_triples, ">", new_alt)
	} else { 
		mutations_snv = NULL
		seq_snvs = c()
	}

	if (nrow(mutations_indel) > 0) {
		# filtering indels in unsequenceable regions of the genome
		indel_ranges = GenomicRanges::GRanges(mutations_indel$chr,
				IRanges::IRanges(start = mutations_indel$pos1, 
								end = mutations_indel$pos2), strand="*")
		indel_seqs = as.character(BSgenome::getSeq(this_genome, indel_ranges))
		seq_indels = grep("N", indel_seqs, invert = T)
		mutations_indel = mutations_indel[seq_indels,]
		indel_seqs = indel_seqs[seq_indels]
		
		mutations_indel$tag = "indel>X"
	} else {
		mutations_indel = NULL
		seq_indels = c()
	}

	# Filtering message
	n_muts_removed = (nrow(mutations) - length(seq_snvs) - length(seq_indels))
	cat("Removing ", n_muts_removed, " invalid SNVs & indels\n\n")
	
	rbind(mutations_snv, mutations_indel)
}


#' This function filters hypermutated samples and returns the formatted mutations with the appropriate trinucleotide context
#'
#' @param filter_hyper_MB The number of mutations per megabase for which a sample is considered hypermutated.
#' Hypermutated samples will be removed in further analyses.
#' @param mutations A data frame with the following columns: chr, pos1, pos2, ref, alt, patient
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{pos1}{the start position of the mutation in base 1 coordinates}
#'     \item{pos2}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#'     \item{patient}{the patient identifier as a string}
#' }
#' 
#' @param this_genome The reference genome object of BSgenome
#'
#' @return a data frame called mutations which has been formatted with an extra column for trinucleotide context
#' @export
#'
#' @examples
#' \donttest{
#' data(cll_mutations)
#' this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#' formatted_mutations = format_muts(cll_mutations[1:10,], 
#' 		filter_hyper_MB = 30, this_genome = this_genome)
#' }
format_muts = function(mutations, this_genome, filter_hyper_MB = NA) {
	
	# remove hypermutated samples, according to muts/megabase rate defined
	if (!is.na(filter_hyper_MB) & filter_hyper_MB > 0) {
		total_muts_filter = 3000*filter_hyper_MB
		sample_mut_count = table(mutations$patient)
		
		hyper_tab = sample_mut_count[sample_mut_count > total_muts_filter]
		spl_rm = names(hyper_tab)
		no_mut_rm = sum(hyper_tab)
		
		cat(length(spl_rm), "remove hypermut, n=", no_mut_rm, ", ", round(100*no_mut_rm/nrow(mutations)), "%\n")
		cat("hypermuted samples: ", spl_rm, "\n\n")
		
		mutations = mutations[!mutations$patient %in% spl_rm,, drop=F]
	}
	if(nrow(mutations) == 0) {
		stop("No mutations left after filtering hypermutators")
	}
#	
#	# keep only relevant chrs, make sure CHR is present in address
#	mutations$chr = toupper(gsub("chr", "", mutations$chr, ignore.case = TRUE))
#	mutations$chr = paste0("chr", mutations$chr)
#	mutations = mutations[mutations$chr %in% paste0("chr", c(1:22, "Y", "X", "M")),]
#	
	mutations$pos1 = as.numeric(mutations$pos1)
	mutations$pos2 = as.numeric(mutations$pos2)
#	
#	if(nrow(mutations) == 0) { 
#		stop("No mutations left after filtering")
#	}
	
	# reverse start/end coordinates of deletions
	rev_coords = which(mutations$pos2 - mutations$pos1 < 0)
	cat("reversing", length(rev_coords), "positions\n")
	if (length(rev_coords)>0) {
		pos1 = mutations[rev_coords, "pos1"]
		pos2 = mutations[rev_coords, "pos2"]
		mutations[rev_coords, "pos1"] = pos2
		mutations[rev_coords, "pos2"] = pos1
		rm(pos1, pos2)
	}
	
	mutations = .get_3n_context_of_mutations(mutations, this_genome)
	mutations
}
