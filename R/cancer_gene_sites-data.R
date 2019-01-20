#' post-translational modification sites found in cancer genes
#'
#' @name cancer_gene_sites
#'
#' @docType data
#'
#' @usage data(cancer_gene_sites)
#'
#' @format A data frame containing the following columns: chr, start, end, id
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the site in base 0 coordinates (BED format)}
#'     \item{end}{the end position of the site in base 0 coordinates (BED format)}
#'     \item{id}{the site identifier - each site should contain only 1 segment and a unique id. If ids are duplicated, each segment of the site will be treated as an individual site. Sites can be coding or noncoding such as phosphosites of protein coding genes in genomic coordinates or transcription factor binding sites of active enhancers.}
#'}
#'
#' @keywords datasets
#'
#' @references Wadi, Lina, et al. "Candidate cancer driver mutations in super-enhancers and long-range chromatin interaction networks." bioRxiv (2017): 236802.
#' (\href{https://www.biorxiv.org/content/early/2017/12/19/236802}{bioRxiv})
#'
#' @source \href{https://www.biorxiv.org/content/early/2017/12/19/236802}{bioRxiv}
#'
#' @examples
#' data(cll_mutations)
#' data(cancer_genes)
#' data(cancer_gene_sites)
#' ActiveDriverWGS(mutations = cll_mutations, elements = cancer_genes, sites = cancer_gene_sites)
NULL
