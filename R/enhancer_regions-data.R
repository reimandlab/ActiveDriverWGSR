#' enhancer regions
#'
#' enhancer regions as referenced in PCAWG
#' 
#' @name enhancer_regions
#' 
#' @docType data
#'
#' @usage data(enhancer_regions)
#'
#' @format A data frame containing the following columns: chr, start, end, id
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
#' @keywords datasets
#'
#' @references Sabarinathan, Radhakrishnan, et al. "The whole-genome panorama of cancer drivers." bioRxiv (2017): 190330.
#' (\href{https://www.biorxiv.org/content/early/2017/09/20/190330}{bioRxiv})
#'
#' @source \href{https://www.biorxiv.org/content/early/2017/09/20/190330}
#'
#' @examples
#' data(pancancer_mutations)
#' data(enhancer_regions)
#' ActiveDriverWGS(mutations = pancancer_mutations, elements = enhancer_regions)
NULL