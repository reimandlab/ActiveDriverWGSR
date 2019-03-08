

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

  # Progress Update
  if (i %% 100 == 0) cat(i, " ")

  # Retrieving coords
  n_regions =  as.numeric(coords[i, "V10"])
  lens = as.numeric(strsplit(coords[i, "V11"], ',')[[1]])
  starts = as.numeric(strsplit(coords[i, "V12"], ',')[[1]]) + coords[i, "V2"]
  ends = starts+lens

  # Checking File Format
  if (length(lens) != n_regions | length(starts) != n_regions) stop("Incorrect BED12 Format: Incorrect blockCounts")
  if (any(starts < coords[i, "V2"]) | any(ends > coords[i, "V3"])) stop("Incorrect BED12 Format: Incorrect blockSizes or blockStarts")

  # Returning Element
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
#' elements = prepare_elements_from_BED12(system.file("extdata",
#' "chr17.coding_regions.bed",
#' package = "ActiveDriverWGS",
#' mustWork = TRUE))
prepare_elements_from_BED12 = function(fname) {

  # Legal Chromosomes
  legal_chr = paste0("chr", c(1:22, "M", "X", "Y"))

  # Reading File
  input = utils::read.delim(fname, stringsAsFactors=F, header=F)
  colnames(input) = paste0("V", 1:ncol(input))

  # Checking File Format
  if (ncol(input) != 12) stop("Incorrect BED12 Format: 12 Columns in BED12 files")
  if (!is.numeric(input$V2) | !is.numeric(input$V3)) stop("Incorrect BED12 Format: Incorrect coordinate format")
  if (!any(input$V1 %in% legal_chr)) stop("Incorrect BED12 Format: Chromosomes must be autosomal, sex or mitochondrial")
  if (!is.character(input$V4)) stop("Incorrect BED12 Format: IDs must be a character string")

  # Processing File
  cat("\n", nrow(input), " Rows :: Processing row ")
  coords = do.call(rbind, lapply(1:nrow(input), .split_coord_fragments_in_BED, input))
  coords = data.frame(coords, stringsAsFactors=F)
  coords$starts = as.numeric(coords$starts)
  coords$ends = as.numeric(coords$ends)
  # coords$frag_id = as.numeric(coords$frag_id) to exclude frag_id as a necessary column
  dim1 = nrow(coords)

  # Filtering Regions for Those Identified on a Single Chromosome (Mitochondrial, Autosomal and Sex Chromosomes)
  # & Repeated Elements
  chr = by(coords$chr, coords$id, function(x) unique(as.character(x)), simplify=F)
  exclude_XY = names(which(sapply(chr, length)==2)) # Multiple chromosomes for 1 ID
  coords = coords[!(coords$id %in% exclude_XY) & coords$chr %in% legal_chr,]
  # remove repeated elements
  coords = unique(coords)
  dim2 = nrow(coords)
  cat("\n Preparing Elements Complete \n RM", dim1-dim2, "lines\n")
  coords
}
