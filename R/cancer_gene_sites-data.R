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
#' @references Wadi, Lina, et al. "ActiveDriverDB: human disease mutations and genome variation in post-translational modification sites of proteins." Nucleic Acids Res. (2018): Jan 4;46(D1):D901-D910.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/29126202/}{PubMed})
#'
#' @source \href{https://pubmed.ncbi.nlm.nih.gov/29126202/}{PubMed}
#'
#' @examples
#' data(cancer_gene_sites)
#' \donttest{
#' data(cll_mutations)
#' data(cancer_genes)
#' ActiveDriverWGS(mutations = cll_mutations, elements = cancer_genes, sites = cancer_gene_sites)
#' }
NULL
