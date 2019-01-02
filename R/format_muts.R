require(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg19)

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
#'
#' @return A data frame consisting of the same columns as the original mutations data frame and sorted
#' by SNVs and Indels with an additional column \code{tag} which indicates the trinucleotide context of
#' the mutation
#' @export
#'
#' @examples
#' get_3n_context_of_mutations(mutations)
get_3n_context_of_mutations = function(mutations) {
  
  legal_dna = c("A", "C", "G", "T")
  
  mutations_snv = mutations[mutations$ref %in% legal_dna & mutations$alt %in% legal_dna,]
  mutations_mnv = mutations[!(mutations$ref %in% legal_dna & mutations$alt %in% legal_dna),]
  
  # snvs can have flanks
  flank_ranges = GRanges(mutations_snv$chr, 
                         IRanges(start=mutations_snv$pos1-1, end=mutations_snv$pos2+1), strand="*")
  triples = as.character(getSeq(Hsapiens, flank_ranges))
  
  # complement trinucleotide where necessary to force into one signature space
  new_triples = triples
  new_alt = mutations_snv$alt
  which_to_complement = which(mutations_snv$ref %in% c("G", "A"))
  
  new_triples[which_to_complement] = as.character(complement(DNAStringSet(new_triples[which_to_complement])))
  new_alt[which_to_complement] = as.character(complement(DNAStringSet(new_alt[which_to_complement])))
  mutations_snv$tag = paste0(new_triples, ">", new_alt)
  #	ref2 = substr(triples, 2,2)
  #	mutations_snv[ref2!=mutations_snv$ref, "tag"] = "dubref"
  #	mutations_snv = mutations_snv[mutations_snv$tag!="dubref",]
  
  if (nrow(mutations_mnv)==0) {
    return(mutations_snv)
  }
  
  mutations_mnv$tag = "indel>X"
  rbind(mutations_snv, mutations_mnv)
}


#' This function filters hypermutated samples and returns the formatted mutations with the appropriate trinucleotide context
#'
#' @param filter_hyper_MB 
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
#' @return a data frame called mutations which has been formatted with an extra column for trinucleotide context
#' @export
#'
#' @examples
format_muts = function(mutations, filter_hyper_MB=NA) {

  # remove hypermutated samples, according to muts/megabase rate defined
  if (!is.na(filter_hyper_MB) & filter_hyper_MB>0) {
    total_muts_filter = 3000*filter_hyper_MB
    sample_mut_count = table(mutations$patient)
    hyper_tab = sample_mut_count[sample_mut_count>total_muts_filter]
    spl_rm = names(hyper_tab)
    no_mut_rm = sum(hyper_tab)
    cat(length(spl_rm), "remove hypermut, n=", no_mut_rm, ", ", round(100*no_mut_rm/nrow(mutations)), "%\n")
    cat("hypermuted samples: ", spl_rm, "\n\n")
    mutations = mutations[!mutations$patient %in% spl_rm,, drop=F]
  }
  
  # keep only relevant chrs, make sure CHR is present in address
  mutations$chr = toupper(gsub("chr", "", mutations$chr, ignore.case=TRUE))
  mutations$chr = paste0("chr", mutations$chr)
  mutations = mutations[mutations$chr %in% paste0("chr", c(1:22, "Y", "X", "M")),]	
  
  mutations$pos1 = as.numeric(mutations$pos1)
  mutations$pos2 = as.numeric(mutations$pos2)
  
  # reverse start/end coordinates of deletions
  rev_coords = which(mutations$pos2-mutations$pos1<0)
  cat("reversing", length(rev_coords), "positions\n")
  if (length(rev_coords)>0) {
    pos1 = mutations[rev_coords, "pos1"]
    pos2 = mutations[rev_coords, "pos2"]
    mutations[rev_coords, "pos1"] = pos2
    mutations[rev_coords, "pos2"] = pos1
    rm(pos1, pos2)
  }
  mutations = get_3n_context_of_mutations(mutations)
  mutations
}