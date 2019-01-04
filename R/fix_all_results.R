
#' Title
#'
#' @param all_results 
#'
#' @return
#' @export
#'
#' @examples
fix_all_results = function(all_results) {
  
  resi = data.frame(all_results, stringsAsFactors=F)
  resi$pp_site = as.numeric(resi$pp_site)
  resi$pp_element = as.numeric(resi$pp_element)
  resi$element_muts_obs = as.numeric(resi$element_muts_obs)
  resi$element_muts_exp = as.numeric(resi$element_muts_exp)
  resi$element_enriched = as.logical(gsub("\\s+", "", resi$element_enriched))
  resi$site_muts_obs = as.numeric(resi$site_muts_obs)
  resi$site_muts_exp = as.numeric(resi$site_muts_exp)
  resi$site_enriched = as.logical(gsub("\\s+", "", resi$site_enriched))
  
  resi[!is.na(resi$pp_element) & resi$pp_element==0,"pp_element"] = 1e-300
  resi[!is.na(resi$pp_site) & resi$pp_site==0,"pp_site"] = 1e-300
  
  resi_tag = resi[,"id"]
  resi = resi[!duplicated(resi_tag),]
  resi
}

#' Title
#'
#' @param all_res 
#'
#' @return
#' @export
#'
#' @examples
get_signf_results = function(all_res) {
  this_results = all_res
  if (nrow(this_results)==0) {
    return(NULL)
  }
  # this is FDR treating element-level NAs as 1s
  this_results$fdr_element = stats::p.adjust(this_results$pp_element, method="fdr", n=nrow(this_results))
  
  this_results = this_results[order(this_results$fdr_element),]
  
  filtered_results = this_results[!is.na(this_results$fdr_element) & this_results$fdr_element<0.05,]
  unsignf_results = this_results[is.na(this_results$fdr_element) | this_results$fdr_element>=0.05,]
  
  if (nrow(filtered_results) + nrow(unsignf_results) != nrow(this_results)) {
    stop("Error: Something unexpected happened when formatting results")
  }
  
  if (nrow(filtered_results)!=0) {
    # site-level FDR perform only on elements with pre-selection of FDR<0.05
    filtered_results$fdr_site = stats::p.adjust(filtered_results$pp_site, method="fdr", n=nrow(filtered_results))
    filtered_results$has_site_mutations = !is.na(filtered_results$fdr_site) & filtered_results$fdr_site<0.05
    filtered_results$has_site_mutations = c("","V")[1+c(filtered_results$has_site_mutations)]
  }
  
  if (nrow(unsignf_results) !=0) {
    unsignf_results$fdr_site = NA
    unsignf_results$has_site_mutations = ""
  }
  
  final_results = rbind(filtered_results, unsignf_results)
  final_results$pp_element = replace(final_results$pp_element, which(is.na(final_results$pp_element)), 1)
  final_results$pp_site = replace(final_results$pp_site, which(is.na(final_results$pp_site)), 1)
  final_results$fdr_element = replace(final_results$fdr_element, which(is.na(final_results$fdr_element)), 1)
  final_results$fdr_site = replace(final_results$fdr_site, which(is.na(final_results$fdr_site)), 1)
  rownames(final_results) = NULL
  final_results
}