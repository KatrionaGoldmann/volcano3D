#' Extract a subset population
#'
#' Subsets data according to the significance groups. 
#' @param polar A polar object including expression data from groups of
#' interest. Created by \code{\link{polar_coords}}.
#' @param significance Which significance factors to subset to. If NULL 
#' levels(syn_polar@polar$sig)[1] is selected. 
#' @param output What object to return. Options are "pvals", "padj", "data",
#'   "df" for subset dataframes, or "polar" to subset the entire 'volc3d' class
#'   object.
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
#' syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
#'                           data = t(syn_example_rld))
#'
#' subset <- significance_subset(syn_polar, "L+", "df")

significance_subset <- function(polar,
                                significance = NULL, 
                                output = "pvalues"){
  
  if(is.null(significance)) significance <- levels(polar@df[[1]]$lab)[1]
  if(! all(significance %in% levels(polar@df[[1]]$lab))){
    stop("`significance` must be in levels(polar@df[[1]]$lab)")
  }
  if(!is(polar, "volc3d")) stop("polar must be a 'volc3d' class object")
  if(! output %in% c("pvals", "padj", "data", "df", "polar")){
    stop("`output` must be one of 'pvals', 'padj', 'data', 'df', 'polar'")
  }
  
  rows <- polar@df[[1]]$lab %in% significance
  polar@pvals <- polar@pvals[rows, ]
  polar@padj <- polar@padj[rows, ]
  polar@data <- polar@data[, rows]
  polar@df[[1]] <- polar@df[[1]][rows, ]
  polar@df[[2]] <- polar@df[[2]][rows, ]
  
  if(output == "polar"){
    return(polar)
  } else {
    return(slot(polar, output))
  }
}

