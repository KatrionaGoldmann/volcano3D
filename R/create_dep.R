setClassUnion("character_or_NULL", c("character", "NULL"))
setClassUnion("df_or_NULL", c("data.frame", "matrix", "NULL"))

#' An S4 class to represent differential expression and pvalues of samples/
#' probes between three groups.
#'
#' @slot pvalues A data frame containing the p-values, adjusted p-values,
#'   and log2(fold changes) for all three comparisons between groups in the 
#'   contrast factor, as well as optional multi-group tests.
#' @slot sampledata Sample data with column ID and contrast
#' @slot contrast The column name in `sampledata`` which contains the
#'   three-group contrast factor used for comparisons.
#' @slot multi_group_test An optional column name prefix for statistical tests 
#' between all three groups
#' @slot expression An optional data frame or matrix containing the
#'   expression data
setClass("dep", slots = list(pvalues = "data.frame",
                             sampledata = "data.frame",
                             contrast = "character",
                             multi_group_test = "character_or_NULL",
                             expression = "df_or_NULL"))

#' Create differential expression pvalues (dep) object
#'
#' This function creates a dep object from a pvalues data frame.
#' @param sampledata A data frame containing the sample information.
#' This must contain: an ID column containing the sample IDs which can be 
#' matched to the expression data and a 
#' contrast column containing the three-level factor used for contrasts.
#' @param contrast The column name in `sampledata` which contains the 
#' three-level factor used for contrast.
#' @param groups The groups to be compared (if NULL this defaults
#' to \code{levels(sampledata[, 'contrasts'])}).
#' @param pvalues A data frame containing three `p_col_suffix` and three
#' `fc_col_suffix` columns: one for each comparison between groups.  Similarly 
#' it can also contain: three optional `padj_col_suffix` columns (if NULL 
#' adjusted p values are calculated using `padjust_method``); and the 'p', 
#' 'padj and 'fc' columns for a three-way test, such as ANOVA or likelihood 
#' ratio test, defined by `multi_group_prefix`.
#' @param expression Optional data frame containing expression data for
#' downstream analysis and visualisation. The rows must contain probes which
#' match the rows in pvalues and the columns must contain samples which match
#' \code{sampledata$ID}.
#' @param p_col_suffix The suffix word to define columns containing p values
#' (default = 'pvalues').
#' @param padj_col_suffix The suffix word to define columns containing adjusted
#' p values (default = 'padj'). If NULL these will be calculated using
#' padjust_method.
#' @param fc_col_suffix The suffix word to define columns containing adjusted
#' p-values (default = 'logFC').
#' @param padjust_method The method used to calculate adjusted p values if 
#' padj_col_suffix is NULL (default = 'BH'). See 
#' \code{\link[stats]{p.adjust}}.
#' @param multi_group_prefix Optional column prefix for statistics (p, padj, 
#' and fold change) across all three groups (typically ANOVA or likelihood 
#' ratio tests). default = NULL.
#' @keywords pvalue
#' @return Returns an S4 dep object containing:
#' \itemize{
#'   \item{'pvalues'} A data frame containing the p-values, adjusted p-values,
#'   and log2(fold changes) for all three comparisons between groups in the 
#'   contrast factor, as well as optional multi-group tests.
#'   \item{'sampledata'} Sample data with column ID and contrast
#'   \item{'contrast'} The column name in `sampledata`` which contains the
#'   three-group contrast factor used for comparisons.
#'   \item{'multi_group_test'} Column name prefix for statistical tests between
#'   all three groups
#'   \item{'expression'} An optional data frame or matrix containing the
#'   expression data
#' }
#' @importFrom stats p.adjust setNames
#' @examples
#' library(volcano3Ddata)
#' data(syn_data)
#' syn_p_obj <- create_dep(sampledata = syn_metadata, 
#'                     contrast = "Pathotype", 
#'                     pvalues = syn_pvalues,
#'                     p_col_suffix ="pvalue", 
#'                     fc_col_suffix = "log2FoldChange",
#'                     multi_group_prefix = "LRT", 
#'                     expression = syn_rld)
#' @references
#'  Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and 
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @export

create_dep <- function(sampledata,
                       contrast,
                       groups = NULL,
                       pvalues,
                       expression = NULL,
                       p_col_suffix = "pvalues",
                       padj_col_suffix = "padj",
                       fc_col_suffix = "logFC",
                       padjust_method = "BH",
                       multi_group_prefix = NULL){
  
  
  # Check for errors
  if(! class(sampledata) %in% c("data.frame")) {
    stop("sampledata must be a data frame")
  }
  if(! contrast %in% colnames(sampledata)) {
    stop("contrast is not a column in sampledata")
  }
  if(! "ID" %in% colnames(sampledata)) {
    stop("There is no ID column in metadata")
  }
  if(! is.null(expression)){
    if(! identical(rownames(expression), rownames(pvalues))){
      stop('The expression row names must be identical to the pvalues row 
           names')
    }
    if(! identical(colnames(expression), as.character(sampledata$ID))) {
      stop('The expression column names must be identical to the sampledata$ID')
    }
  }
  
  # Ensure groups and contrast column are compatible
  sampledata[, contrast] <- droplevels(factor(sampledata[, contrast]))
  if(length(levels(sampledata[, contrast])) != 3) {
    stop("There number of factors in the comparison column does not equal 3")
  }
  if(is.null(groups)) {groups <- levels(sampledata[, contrast])}
  if(length(groups) != 3) stop("There number of groups does not equal 3")
  if(! is.null(groups)) {if(! all(groups %in% levels(sampledata[, contrast]))) {
    stop('Make sure all groups are in sampledata[, contrast]')
  }
  }
  
  comparisons <- c(paste(groups[1], groups[2], sep="-"),
                   paste(groups[2], groups[3], sep="-"),
                   paste(groups[3], groups[1], sep="-"))
  
  
  # Check column names of correct format exist in pvalues data frame
  for(col_suffix in c(p_col_suffix, fc_col_suffix, padj_col_suffix)){
    comparitiveCols <- paste(c(comparisons, multi_group_prefix), col_suffix)
    notFinding <- c()
    if(! all(comparitiveCols %in% colnames(pvalues))) {
      notFinding <- comparitiveCols[! comparitiveCols %in% colnames(pvalues)]
      notFinding <- notFinding[! is.na(notFinding)]
      
      # check if ordering of groups in column names is the wrong way round
      check <- strsplit(gsub(" ", "", gsub(col_suffix, "", notFinding)), "-")
      
      for (order_check in check){
        og <- paste0(order_check[1], "-", order_check[2], " ", col_suffix)
        reverse <- paste0(order_check[2], "-", order_check[1], " ", col_suffix)
        if(reverse %in% colnames(pvalues)){
          colnames(pvalues)[colnames(pvalues) == reverse] <- og
          
          # Need to reverse order for fold change
          if(col_suffix == fc_col_suffix) {
            pvalues[, og] <- -1*pvalues[, og]
          }
          warning(paste(og, "was not found in colnames(pvalues), but", reverse, 
                        "was - the column name has now been reversed"))
          notFinding <- notFinding[notFinding != og]
        }
      }
    }
    
    if(length(notFinding) > 0){  
      if(length(notFinding) > 1) {
        notFinding <- paste0("'", 
                             paste0(notFinding[1:(length(notFinding)-1)],
                                    collapse="', '"),
                             "' or '", 
                             notFinding[length(notFinding)], "'")
      }
      stop(paste('Cannot find', paste0(notFinding, collapse=", "),
                 'in colnames(pvalues)'))
    }
  }
  
  # If adjusted p is not calculated, calculate
  if(is.null(padj_col_suffix)) {
    for(comp in c(comparisons, multi_group_prefix)){
      pvalues$new <- p.adjust(pvalues[, paste(comp, p_col_suffix)],
                              method = padjust_method)
      colnames(pvalues)[colnames(pvalues) == "new"] <- paste(comp, "padj")
    }
    padj_col_suffix <- "padj"
  }
  
  pOutput <- pvalues[, sort(paste(c(comparisons, multi_group_prefix),
                                 rep(c(p_col_suffix, fc_col_suffix,
                                       padj_col_suffix),
                                     each = length(c(comparisons,
                                                     multi_group_prefix)))))]
  
  colnames(pOutput) <- gsub(p_col_suffix, "pvalue", colnames(pOutput))
  colnames(pOutput) <- gsub(fc_col_suffix, "logFC", colnames(pOutput))
  colnames(pOutput) <- gsub(padj_col_suffix, "padj", colnames(pOutput))
  
  methods::new("dep",
      pvalues = pOutput,
      sampledata = sampledata,
      contrast = contrast,
      multi_group_test = multi_group_prefix,
      expression = expression
  )
  
}



