

#' Splits a BED12 file into separate regions
#'
#' @param i The i-th row of the coords data frame which needs to be split into separate elements
#' @param coords The coords data frame which is the imported BED12 file
#'
#' @return A data frame containing the following columns for a given BED12 identifier
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
.split_coord_fragments_in_BED = function(i, coords) {
  if (i %% 100 == 0) cat(i, " ")
  lens = as.numeric(strsplit(coords[i, "V11"], ',')[[1]])
  starts = as.numeric(strsplit(coords[i, "V12"], ',')[[1]]) + coords[i, "V2"]
  ends = starts+lens
  
  dfr = cbind(chr=coords[i,"V1"], starts, ends, id=coords[i,"V4"]) # , frag_id=1:length(lens)) to exclude frag_id as a necessary column
  dfr	
}


#' Prepares element coords from a BED12 file
#'
#' @param fname The file name of a BED12 file containing the desired elements. For further documentation
#' on the BED12 format, refer to the UCSC website.
#'
#' @return A data frame containing the following columns to be used as the input element coords
#' to ActiveDriverWGS
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the element in base 0 coordinates (BED format)}
#'     \item{end}{the end position of the element in base 0 coordinates (BED format)}
#'     \item{id}{the element identifier - if the element contains multiple segments such as exons,
#'     each segment should be a separate row with the segment coordinates
#'     and the element identifier as id. Elements can be coding or noncoding such as exons of protein
#'     coding genes or active enhancers.}
#' }
#' @export
#'
#' @examples
#' # Don't run
#' elements = prepare_elements_from_BED12("gc19_pc.cds.bed")
prepare_elements_from_BED12 = function(fname) {
  input = read.delim(fname, stringsAsFactors=F, header=F)
  colnames(input) = paste0("V", 1:ncol(input))
  cat("\n", nrow(input), "::  ")
  coords = do.call(rbind, lapply(1:nrow(input), .split_coord_fragments_in_BED, input))
  coords = data.frame(coords, stringsAsFactors=F)
  coords$starts = as.numeric(coords$starts)
  coords$ends = as.numeric(coords$ends)
  coords$frag_id = as.numeric(coords$frag_id)
  dim1 = nrow(coords)
  
  chr = by(coords$chr, coords$id, function(x) unique(as.character(x)), simplify=F)
  exclude_XY = names(which(sapply(chr, length)==2))
  legal_chr = paste0("chr", c(1:22, "M", "X", "Y"))
  coords = coords[!(coords$id %in% exclude_XY) & coords$chr %in% legal_chr,]
  # remove repeated elements
  coords = unique(coords)
  dim2 = nrow(coords)
  cat("RM", dim1-dim2, "lines\n")
  coords
}
