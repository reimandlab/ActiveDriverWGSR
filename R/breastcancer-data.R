#' breast cancer mutations
#'
#' breast cancer whole genome simple somatic mutations from Alexandrov et, 2013
#' 
#' @name breastcancer_mutations
#' 
#' @docType data
#'
#' @usage data(breastcancer_mutations)
#'
#' @format A data frame containing the following columns: chr, pos1, pos2, ref, alt, patient.
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{pos1}{the start position of the mutation in base 1 coordinates}
#'     \item{pos2}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#'     \item{patient}{the patient identifier as a string}
#'     \item{cc}{cancer type of the sample}
#'
#' @keywords datasets
#'
#' @references Alexandrov, Ludmil B., et al. "Signatures of mutational processes in human cancer." Nature 500.7463 (2013): 415.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/23945592}{PubMed})
#'
#' @source \href{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl}{FTP Server}
#'
#' @examples
#' data(breastcancer_mutations)
#' data(cds_coords)
#' ActiveDriverWGS(mutations = breastcancer_mutations, elements = cds_coords)
NULL