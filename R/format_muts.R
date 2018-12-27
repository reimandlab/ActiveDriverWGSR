library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

#' Title
#'
#' @param mutations 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param fname 
#' @param cols_we_need 
#' @param sample_blacklist 
#' @param sample_name_col 
#' @param read_nrow 
#' @param filter_hyper_MB 
#'
#' @return
#' @export
#'
#' @examples
format_muts = function(fname, cols_we_need, sample_blacklist, sample_name_col="file.name", 
                      read_nrow=100000000, filter_hyper_MB=NA) {
  maf = read.delim(fname, stringsAsFactors=F, nrow=read_nrow)
  if (nrow(maf)==read_nrow) {
    stop(paste("warning:", fname, "was not read completely\n"))		
  }
  
  if (!all(cols_we_need %in% colnames(maf))) {
    stop(paste("error:", fname, "has wrong colnames\n"))		
  }
  
  # remove blacklisted samples
  which_black = which(maf[, sample_name_col] %in% sample_blacklist)
  cat("blacklisted", length(which_black), "\n", sample_blacklist, "\n")
  maf = maf[!maf[, sample_name_col] %in% sample_blacklist,]
  
  # keep only needed columns
  maf = maf[,cols_we_need]
  colnames(maf) = c("chr", "pos1", "pos2", "ref", "alt", "patient", "cc")
  gc()
  
  # remove hypermutated samples, according to muts/megabase rate defined
  if (!is.na(filter_hyper_MB) & filter_hyper_MB>0) {
    total_muts_filter = 3000*filter_hyper_MB
    sample_mut_count = table(maf$patient)
    hyper_tab = sample_mut_count[sample_mut_count>total_muts_filter]
    spl_rm = names(hyper_tab)
    no_mut_rm = sum(hyper_tab)
    cat(length(spl_rm), "remove hypermut, n=", no_mut_rm, ", ", round(100*no_mut_rm/nrow(maf)), "%\n")
    cat("hypermuted samples: ", spl_rm, "\n\n")
    maf = maf[!maf$patient %in% spl_rm,, drop=F]
  }
  
  # keep only relevant chrs, make sure CHR is present in address
  maf$chr = toupper(gsub("chr", "", maf$chr, ignore.case=TRUE))
  maf$chr = paste0("chr", maf$chr)
  maf = maf[maf$chr %in% paste0("chr", c(1:22, "Y", "X", "M")),]	
  
  maf$pos1 = as.numeric(maf$pos1)
  maf$pos2 = as.numeric(maf$pos2)
  
  # reverse start/end coordinates of deletions
  rev_coords = which(maf$pos2-maf$pos1<0)
  cat("reversing", length(rev_coords), "positions\n")
  if (length(rev_coords)>0) {
    pos1 = maf[rev_coords, "pos1"]
    pos2 = maf[rev_coords, "pos2"]
    maf[rev_coords, "pos1"] = pos2
    maf[rev_coords, "pos2"] = pos1
    rm(pos1, pos2)
  }
  maf = get_3n_context_of_mutations(maf)
  maf
}