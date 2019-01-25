#' cancer_genes
#'
#' protein coding genes from gencode v.19, cancer genes adapted from the Cancer Gene Census (November, 2018).
#' Genes affected solely by amplifications, deletions and translations were removed.
#'
#' @name cancer_genes
#'
#' @docType data
#'
#' @usage data(cancer_genes)
#'
#' @format A data frame containing the following columns: chr, start, end, id
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the element in base 0 coordinates (BED format)}
#'     \item{end}{the end position of the element in base 0 coordinates (BED format)}
#'     \item{id}{the element identifier - if the element contains multiple segments such as exons, each segment should be a separate row with the segment coordinates and the element identifier as id. Elements can be coding or noncoding such as exons of protein coding genes or active enhancers.}
#' }
#'
#' @keywords datasets
#'
#' @references Harrow, Jennifer, et al. "GENCODE: the reference human genome annotation for The ENCODE Project." Genome research 22.9 (2012): 1760-1774.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/22955987}{PubMed})
#'
#' @source \href{https://www.gencodegenes.org/human/release_19.html}{GENCODE}
#'
#' @examples
#' data(cll_mutations)
#' data(cancer_genes)
#' ActiveDriverWGS(mutations = cll_mutations, elements = cancer_genes)
NULL
