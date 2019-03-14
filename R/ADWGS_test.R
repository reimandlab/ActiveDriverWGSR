#' Makes mutational signatures
#'
#' @return a dataframe with mutational signatures
.make_mut_signatures = function() {
  nucs = c("A", "C", "G", "T")

  signatures_dfr = data.frame(as.matrix(expand.grid(left=nucs, mid=nucs, right=nucs, alt=nucs)), stringsAsFactors=F)
  signatures_dfr = signatures_dfr[signatures_dfr$mid!=signatures_dfr$alt,]
  signatures_dfr = signatures_dfr[signatures_dfr$mid %in% c("C", "T"),]
  signatures_dfr$tri = paste0(signatures_dfr$left, signatures_dfr$mid, signatures_dfr$right)
  signatures_dfr$signt = paste0(signatures_dfr$tri, ">", signatures_dfr$alt)
  signatures_dfr
}

#' Calculates the number of expected mutations based
#'
#' @param hyp hypothesis to be tested
#' @param select_positions boolean column which indicates which positions are in the element of interest
#' @param dfr a dataframe containing the data to be tested
#' @param colname name of the column which indicates the count of mutations in the positions of interest
#'
#' @return a list of observed mutations (numeric), expected mutations (numeric),
#' observations enriched (boolean) and observations depleted (boolean)
.get_obs_exp = function(hyp, select_positions, dfr, colname) {
  obs_mut = sum(dfr[select_positions, colname])

  exp_probs = hyp$fitted.values[select_positions]
  # simulate from poisson distribution with values from model response
  exp_boot = replicate(1000, sum(stats::rpois(rep(1, length(exp_probs)), exp_probs)))
  exp_mut = stats::median(exp_boot)
  obs_enriched = obs_mut > stats::quantile(exp_boot, 0.95)
  obs_depleted = obs_mut < stats::quantile(exp_boot, 0.05)

  list(obs_mut, exp_mut, obs_enriched, obs_depleted)
}

# @import GenomicRanges
# @import GenomeInfoDb
# @import BSgenome.Hsapiens.UCSC.hg19
# @import IRanges
# @import BSgenome
# @import S4Vectors
# @import plyr

#' ADWGS_test executes the statistical test for ActiveDriverWGS
#'
#' @param id A string used to identify the element of interest. \code{id}
#' corresponds to an element in the id column of the elements file
#' @param gr_element_coords A GenomicRanges object that describes the elements of interest containing the
#' chromosome, start and end coordinates, and an mcols column corresponding to id
#' @param gr_site_coords  A GenomicRanges object that describes the sites of interest which reside
#' in the elements of interest containing the chromosome, start and end coordinates,
#' and an mcols column corresponding to id. Examples of sites include transcription factor binding
#' sites in promoter regions or phosphosites in exons of protein coding genes. An empty GenomicRanges object
#' nullifies the requirement for sites to exist.
#' @param gr_maf A GenomicRanges object that describes the mutations in the dataset containing the chromosome,
#' start and end coordinates, patient id, and trinucleotide context
#' @param win_size An integer indicating the size of the background window in base pairs for which the
#' mutation rate is expected to remain the same. The default is 50000bps.
#' @param element_bias A boolean indicating whether or not indels should be counted by their midpoints
#' or with bias towards the element
#'
#' @return A data frame containing the following columns
#' \describe{
#'     \item{id}{A string identifying the element of interest}
#'     \item{pp_element}{The p-value of the element}
#'     \item{element_muts_obs}{The number of patients with a mutation in the element}
#'     \item{element_muts_exp}{The expected number of patients with a mutation in the element with respect to background}
#'     \item{element_enriched}{A boolean indicating whether the element is enriched in mutations}
#'     \item{pp_site}{The p-value of the site}
#'     \item{site_muts_obs}{The number of patients with a mutation in the site}
#'     \item{site_muts_exp}{The expected number of patients with a mutation in the site with respect to element}
#'     \item{site_enriched}{A boolean indicating whether the site is enriched in mutations}
#'     \item{result_number}{A numeric indicator denoting the order in which the results were calculated}
#'     \item{fdr_element}{The FDR corrected p-value of the element}
#'     \item{fdr_site}{The FDR corrected p-value of the site}
#'     \item{has_site_mutations}{A V indicates the presence of site mutations}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(GenomicRanges)
#'
#' # Regions
#' data(cancer_genes)
#' gr_element_coords = GRanges(seqnames = cancer_genes$chr,
#' IRanges(start = cancer_genes$start, end = cancer_genes$end),
#' mcols = cancer_genes$id)
#'
#' # Sites (NULL)
#' gr_site_coords = GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
#'
#' # Mutations
#' data(cll_mutations)
#' cll_mutations = format_muts(cll_mutations)
#'
#' gr_maf = GRanges(cll_mutations$chr,
#' IRanges(cll_mutations$pos1, cll_mutations$pos2),
#' mcols=cll_mutations[,c("patient", "tag")])
#'
#' # ADWGS_test
#' id = "ATM"
#' result = ADWGS_test(id, gr_element_coords, gr_site_coords, gr_maf, 50000)
#'}
ADWGS_test = function(id, gr_element_coords, gr_site_coords, gr_maf, win_size, element_bias = T) {

  cat(".")
  null_res = data.frame(id,
                        pp_element=NA, element_muts_obs=NA, element_muts_exp=NA, element_enriched=NA,
                        pp_site=NA, site_muts_obs=NA, site_muts_exp=NA, site_enriched=NA,
                        stringsAsFactors=F)

  gr_element = gr_element_coords[GenomicRanges::mcols(gr_element_coords)[,1]==id]

  # background is plus/minus window and takes care of cases on chromosome borders
  background_chr = as.character(GenomicRanges::seqnames(gr_element))[1]
  background_start = min(GenomicRanges::start(gr_element))-win_size-2
  background_start = max(background_start, 2)
  background_end = max(GenomicRanges::end(gr_element))+win_size+1
  background_end = min(background_end, GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[background_chr])
  gr_background_plus_element = GenomicRanges::GRanges(background_chr,
                                                      IRanges::IRanges(background_start, background_end))

  # get sequence trinucleotide
  this_seq = strsplit(as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, gr_background_plus_element)), '')[[1]]
  left_nucleotide = this_seq[1:(length(this_seq)-2)]
  mid_nucleotide = this_seq[2:(length(this_seq)-1)]
  right_nucleotide = this_seq[3:(length(this_seq)+0)]

  # reverse strand where needed
  reverser = c("A" = "T", "G" = "C", "T" = "A", "C" = "G", "N" = "N")
  where_reverse = mid_nucleotide %in% c("A", "G")

  # tri nucleotide takes left and right neighbour
  mid_nucleotide[where_reverse] = reverser[mid_nucleotide[where_reverse]]
  left_nucleotide[where_reverse] = reverser[left_nucleotide[where_reverse]]
  right_nucleotide[where_reverse] = reverser[right_nucleotide[where_reverse]]
  tri_nucleotide = paste0(left_nucleotide, mid_nucleotide, right_nucleotide)
  rm(mid_nucleotide, left_nucleotide, right_nucleotide, this_seq, where_reverse)

  # background is plus/minus window size around each segment of element
  # in some cases intermediate non-segments (introns) extend well beyond window size. remove these remote coordinates from background
  exon_bg_posi = unique(unlist(lapply(1:length(gr_element),
                                      function(i) (GenomicRanges::start(gr_element)[i]-win_size):(GenomicRanges::end(gr_element)[i]+win_size))))

  # data frame with position and quad-nucleotide context
  mut_signt = .make_mut_signatures()
  positions = (GenomicRanges::start(gr_background_plus_element)+1):(GenomicRanges::end(gr_background_plus_element)-1)
  dfr = data.frame(pos=positions, nucl=tri_nucleotide, stringsAsFactors=F)
  dfr = dfr[dfr$pos %in% exon_bg_posi,, drop=F]
  dfr = merge(dfr, mut_signt[,c("tri", "signt")], by.x="nucl", by.y="tri")
  dfr = rbind(dfr,
              data.frame(nucl="indel", pos=unique(dfr$pos), signt="indel>X", stringsAsFactors=F))
  dfr = dfr[order(dfr$pos),]
  rm(exon_bg_posi, positions, tri_nucleotide)

  # capture element mutations first and remove repeated mutations per patient
  gr_maf_element = gr_maf[S4Vectors::queryHits(suppressWarnings(GenomicRanges::findOverlaps(gr_maf, gr_element)))]
  gr_maf_element = gr_maf_element[!duplicated(GenomicRanges::mcols(gr_maf_element)[,1])]
  if (length(gr_maf_element)==0) {
    return(null_res)
  }

  # capture background mutations
  gr_background_only = GenomicRanges::setdiff(gr_background_plus_element, gr_element)
  gr_maf_background_only = gr_maf[S4Vectors::queryHits(suppressWarnings(GenomicRanges::findOverlaps(gr_maf, gr_background_only, type = "within")))]

  # count mutations by position and quad-nucl context; indels will be counted by midpoint
  gr_maf_element_plus_background = c(gr_maf_element, gr_maf_background_only)
  # gr_maf_element_plus_background = gr_maf_element_plus_background[!duplicated(data.frame(gr_maf_element_plus_background))]

  # Counting all mutations which overlap the element as element_mutations
  if(element_bias){
    this_gr_maf_element_plus_background = GenomicRanges::granges(gr_maf_element_plus_background)
    element_edge_indels = intersect(S4Vectors::queryHits(suppressWarnings(GenomicRanges::findOverlaps(this_gr_maf_element_plus_background, gr_element))),
                                    S4Vectors::queryHits(suppressWarnings(GenomicRanges::findOverlaps(this_gr_maf_element_plus_background, gr_background_only))))
    if(length(element_edge_indels) > 0)
      this_gr_maf_element_plus_background[element_edge_indels] = do.call(c, lapply(element_edge_indels, function(i) {GenomicRanges::setdiff(this_gr_maf_element_plus_background[i], gr_background_only)}))
    S4Vectors::mcols(this_gr_maf_element_plus_background) = S4Vectors::mcols(gr_maf_element_plus_background)
    gr_maf_element_plus_background = this_gr_maf_element_plus_background
    rm(this_gr_maf_element_plus_background)
  }

  # Creating the muts_per_pos data.frame where there are columns for position, tag and number of mutations for background and element
  mut_midpoint = round((GenomicRanges::start(gr_maf_element_plus_background) + GenomicRanges::end(gr_maf_element_plus_background))/2)
  mut_tag = gr_maf_element_plus_background$mcols.tag
  maf_element_plus_background = data.frame(pos1=mut_midpoint, tag=mut_tag, stringsAsFactors=F)
  muts_per_pos = plyr::ddply(maf_element_plus_background, c("pos1", "tag"), function(x) nrow(x))
  colnames(muts_per_pos)[3] = "n_mut"

  # combine entire region and its mutations by positions and nucleotide content
  dfr = merge(dfr, muts_per_pos, by.x=c("pos", "signt"), by.y=c("pos1", "tag"), all=T)
  dfr[is.na(dfr$n_mut), "n_mut"] = 0

  # label nucleotides that are part of the element
  element_posi = unique(unlist(lapply(1:length(gr_element), function(i) GenomicRanges::start(gr_element)[i]:(GenomicRanges::end(gr_element)[i]))))
  dfr$is_element = dfr$pos %in% element_posi

  # label nucleotides that are part of sites in the element
  gr_element_sites = gr_site_coords[S4Vectors::subjectHits(suppressWarnings(GenomicRanges::findOverlaps(gr_element, gr_site_coords)))]
  site_posi = c()
  if (length(gr_element_sites)>0) {
    site_posi = unique(unlist(lapply(1:length(gr_element_sites),
                                     function(i) GenomicRanges::start(gr_element_sites)[i]:(GenomicRanges::end(gr_element_sites)[i]))))
    site_posi = intersect(site_posi, element_posi)
  }
  dfr$is_site = dfr$pos %in% site_posi

  # simplification of models - discarding nucleotide context with zero mutations throughout
  silent_signt = names(which(by(dfr$n_mut, dfr$signt, sum)==0))
  dfr = dfr[!dfr$signt %in% silent_signt,]
  dfr = dfr[, !colnames(dfr) %in% c("nucl", "pos")]

  # test mutation rate hypothesis
  h0_formula = "n_mut~signt"
  if (length(unique(dfr$signt)) == 1) {
    h0_formula = "n_mut~1"
  }

  h0 = stats::glm(stats::as.formula(h0_formula), data=dfr, family=stats::poisson)
  h1 = stats::update(h0, . ~ . + is_element)
  pp_element = stats::anova(h0, h1, test="Chisq")[2,5]

  element_stats = .get_obs_exp(h0, dfr$is_element, dfr, "n_mut")
  element_muts_obs = element_stats[[1]]
  element_muts_exp = element_stats[[2]]
  element_enriched = element_stats[[3]]
  element_depleted = element_stats[[4]]
  rm(h0)

  # p-values significantly depleted elements should be inverted
  if (element_depleted & !is.na(pp_element) & pp_element<0.5) {
    pp_element = 1 - pp_element
  }

  # if region has sites, test second hypothesis on regions
  pp_site = site_muts_obs = site_muts_exp = site_enriched = site_depleted = NA
  if (sum(dfr$is_site)>0){
    h2 = stats::update(h1, . ~ . + is_site)
    pp_site = stats::anova(h1, h2, test="Chisq")[2,5]

    site_stats = .get_obs_exp(h1, dfr$is_site, dfr, "n_mut")
    site_muts_obs = site_stats[[1]]
    site_muts_exp = site_stats[[2]]
    site_enriched = site_stats[[3]]
    site_depleted = site_stats[[4]]
    rm(h2)

    if (site_depleted & !is.na(pp_site) & pp_site<0.5) {
      pp_site = 1 - pp_site
    }
  }

  rm(h1, dfr)

  data.frame(id,
             pp_element, element_muts_obs, element_muts_exp, element_enriched,
             pp_site, site_muts_obs, site_muts_exp, site_enriched,
             stringsAsFactors=F)
}
