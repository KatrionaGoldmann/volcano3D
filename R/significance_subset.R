#' Extract a subset population
#'
#' Subsets data according to the significance groups. 
#' @param polar A polar object including expression data from groups of
#' interest. Created by \code{\link{polar_coords}}.
#' @param significance Which significance factors to subset to. If NULL 
#' levels(syn_polar@polar$sig)[1] is selected. 
#' @param output What object to return. Options are "pvalues", "expression", 
#' "polar_df" for subset data frames or "polar" for subset polar class object. 
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @importFrom methods slot
#' @export
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(sampledata = syn_example_meta,
#'                           contrast = "Pathotype",
#'                           groups = NULL,
#'                           pvalues = syn_example_p,
#'                           expression = syn_example_rld,
#'                           p_col_suffix = "pvalue",
#'                           padj_col_suffix = "padj",
#'                           non_sig_name = "Not Significant",
#'                           multi_group_prefix = "LRT",
#'                           significance_cutoff = 0.01,
#'                           fc_col_suffix='log2FoldChange',
#'                           fc_cutoff = 0.3)
#'
#' subset <- significance_subset(syn_polar, "Lymphoid+", "polar_df")

significance_subset <- function(polar,
                                significance = NULL, 
                                output = "pvalues"){
  
  if(is.null(significance)) significance <- levels(polar@polar$sig)[1]
  if(! all(significance %in% levels(polar@polar$sig))){
    stop("Significance must be in levels(polar@polar$sig)")
  }
  if(class(polar) != "polar") stop("polar must be a polar class object")
  if(! output %in% c("pvalues", "expression", "polar_df", "polar")){
    stop("return must be one of 'pvalues', 'expression', 'polar_df', 'polar'")
  }
  
  rows <- polar@polar$Name[polar@polar$sig %in% significance]
  polar@pvalues <- polar@pvalues[rows, ]
  polar@expression <- polar@expression[rows, ]
  polar@polar <- polar@polar[rows, ]
  
  if(output == "polar"){
    return(polar)
  } else {
    return(slot(polar, gsub("_df", "", output)))
  }
}

