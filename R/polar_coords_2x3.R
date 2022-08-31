
#' Coordinates for three way polar plot from 2x3 factor analysis
#'
#' This function creates a 'volc3d' object of S4 class for downstream plots
#' containing the p-values from a 2x3 factor analysis, expression data
#' sample data and polar coordinates. For RNA-Seq count data, two functions
#' \code{\link{deseq_2x3}} followed by [deseq_2x3_polar] can be used instead.
#'
#' @param data Dataframe or matrix with variables in columns and samples in rows
#' @param metadata Dataframe of sample information with samples in rows
#' @param outcome Either the name of column in `metadata` containing the binary
#'   outcome data. Or a vector with 2 groups, ideally a factor. If it is not a
#'   factor, this will be coerced to a factor. This must have exactly 2 levels.
#' @param group Either the name of column in `metadata` containing the 3-way
#'   grouping data. Or a vector with 3 groups, ideally a factor. If it is not a
#'   factor, this will be coerced to a factor. This must have exactly 3 levels.
#'   NOTE: if `pvals` is given, the order of the levels in `group` must
#'   correspond to the order of columns in `pvals`.
#' @param pvals Optional matrix or dataframe with p-values in 3 columns. If
#'   `pvals` is not given, it is calculated using the function
#'   \code{\link{calc_stats_2x3}}. The p-values in 3 columns represent the
#'   comparison between the binary outcome with each column for the 3 groups as
#'   specified in `group`.
#' @param padj Matrix or dataframe with adjusted p-values. If not supplied,
#'   defaults to use nominal p-values from `pvals`.
#' @param pcutoff Cut-off for p-value significance
#' @param padj.method Can be `"qvalue"` or any method available in `p.adjust`.
#'   The option `"none"` is a pass-through.
#' @param process Character value specifying colour process for statistical
#'   significant genes: "positive" specifies genes are coloured if fold change
#'   is >0; "negative" for genes with fold change <0 (note that for clarity the
#'   polar position is altered so that genes along each axis have the most
#'   strongly negative fold change values); or "two.sided" which is a compromise
#'   in which positive genes are labelled as before but genes with negative fold
#'   changes and significant p-values have an inverted colour scheme.
#' @param scheme Vector of colours starting with non-significant variables
#' @param labs Optional character vector for labelling groups. Default `NULL`
#'   leads to abbreviated labels based on levels in `outcome` using
#'   [abbreviate()]. A vector of length 3 with custom abbreviated names for the
#'   outcome levels can be supplied. Otherwise a vector length 7 is expected, of
#'   the form "ns", "B+", "B+C+", "C+", "A+C+", "A+", "A+B+", where "ns" means
#'   non-significant and A, B, C refer to levels 1, 2, 3 in `outcome`, and must
#'   be in the correct order.
#' @param ... Optional arguments passed to \code{\link{calc_stats_2x3}}
#' @details
#' This function is designed for manually generating a 'volc3d' class object for
#' visualising a 2x3 way analysis comparing a large number of attributes such as
#' genes. For RNA-Seq data we suggest using [deseq_2x3] and [deseq_2x3_polar]
#' functions in sequence instead.
#' 
#' Scaled polar coordinates are generated using the t-score for each group
#' comparison. Unscaled polar coordinates are generated as difference between
#' means for each group comparison. If p-values are not supplied they are
#' calculated by [calc_stats_2x3] using either t-tests or wilcoxon tests.
#' 
#' The colour scheme is not as straightforward as for the standard polar plot
#' and volcano3D plot since genes (or attributes) can be significantly up or
#' downregulated in the response comparison for each of the 3 groups.
#' `process = "positive"` means that genes are labelled with colours if a gene
#' is significantly upregulated in the response for that group. This uses the
#' primary colours (RGB) so that if a gene is upregulated in both red and blue
#' group it becomes purple etc with secondary colours. If the gene is
#' upregulated in all 3 groups it is labelled black. Non-significant genes are
#' in grey.
#' 
#' With `process = "negative"` genes are coloured when they are significantly
#' downregulated. With `process = "two.sided"` the colour scheme means that both
#' significantly up- and down-regulated genes are coloured with downregulated
#' genes labelled with inverted colours (i.e. cyan is the inverse of red etc).
#' However, significant upregulation in a group takes precedence.
#' 
#' @return Returns an S4 'volc3d' object containing:
#' \itemize{
#'   \item{'df'} A list of 2 dataframes. Each dataframe contains both x,y,z
#'   coordinates as well as polar coordinates r, angle. The first dataframe has
#'   coordinates on scaled data. The 2nd dataframe has unscaled data (e.g. log2
#'   fold change for gene expression). The `type` argument in
#'   \code{\link{volcano3D}}, \code{\link{radial_plotly}} and
#'   \code{\link{radial_ggplot}} corresponds to these dataframes.
#'   \item{'outcome'} The three-group contrast factor used for comparisons,
#'   linked to the `group` column
#'   \item{'data'} Dataframe or matrix containing the expression data
#'   \item{'pvals'} A dataframe containing p-values in 3 columns representing
#'   the binary comparison for the outcome for each of the 3 groups.
#'   \item{'padj'} A dataframe containing p-values adjusted for multiple testing
#'   \item{'pcutoff} Numeric value for cut-off for p-value significance
#'   \item{'scheme'} Character vector with colour scheme for plotting
#'   \item{'labs'} Character vector with labels for colour groups
#' }
#' 
#' @seealso \code{\link{deseq_2x3}}, \code{\link{deseq_2x3_polar}},
#'   \code{\link{calc_stats_2x3}}
#' @export

polar_coords_2x3 <- function(data,
                             metadata = NULL,
                             outcome,
                             group,
                             pvals = NULL, 
                             padj = pvals,
                             pcutoff = 0.05,
                             padj.method = "BH",
                             process = c("two.sided", "positive", "negative"),
                             scheme = c('grey60', 'red', 'gold2', 'green3', 
                                        'cyan', 'blue', 'purple', 'black'),
                             labs = NULL,
                             ...) {
  process <- match.arg(process)
  if (!is.null(metadata)) {
    if (length(outcome) != 1 | !outcome %in% colnames(metadata)) {
      stop(outcome, " is not a column in metadata")
    }
    outcome <- factor(metadata[, outcome])
    if (length(group) != 1 | !group %in% colnames(metadata)) {
      stop(group, " is not a column in metadata")
    }
    group <- factor(metadata[, group])
  }
  if (!all.equal(nrow(data), length(outcome), length(group))) {
    stop("Mismatch between data and metadata")
  }
  if (nlevels(group) != 3) stop("Number of levels in group is not 3")
  if (nlevels(outcome) != 2) stop("Number of levels in outcome is not 2")
  
  res <- calc_stats_2x3(data, outcome, group, pcutoff, padj.method, ...)
  rn <- Reduce(intersect, lapply(res, rownames))
  df1 <- getCols(res, rn, 'statistic')
  df2 <- getCols(res, rn, 'mean.diff')
  if (is.null(pvals)) {
    pvals <- getCols(res, rn, 'pvalue')
    padj <- getCols(res, rn, 'padj')
  }
  SE <- getCols(res, rn, 'stderr')
  colnames(SE) <- paste0("SE_", colnames(SE))
  if (process == "negative") {
    df1 <- -df1
    df2 <- -df2
  }
  df1 <- cbind(df1, SE)
  df2 <- cbind(df2, SE)
  df1 <- polar_xy(df1)
  df2 <- polar_xy(df2)
  outcome <- factor(levels = levels(group))
  
  ptab <- polar_p2x3(df2, pvals, padj, pcutoff, process, scheme, labs)
  df1 <- cbind(df1, ptab)
  df2 <- cbind(df2, ptab)
  
  methods::new("volc3d",
               df = list(scaled = df1, unscaled = df2, type = "polar_coords_2x3"),
               outcome = outcome,
               data = data, pvals = pvals, padj = padj,
               pcutoff = pcutoff, scheme = scheme,
               labs = levels(ptab$lab))
}


#' Calculate p-values for 2x3-way analysis
#' 
#' @param data Dataframe or matrix with variables in columns and samples in rows
#' @param outcome Factor with 2 levels
#' @param group Factor with 3 levels
#' @param pcutoff Cut-off for p-value significance
#' @param padj.method Can be `"qvalue"` or any method available in `p.adjust`.
#'   The option `"none"` is a pass-through.
#' @param test Character value specifying the statistical test between the 2
#'   level response outcome. Current options are "t.test" or "wilcoxon".
#' @param exact Logical for whether to use an exact test (Wilcoxon test only)
#' @importFrom matrixTests col_t_welch col_wilcoxon_twosample
#' @export
#' 
calc_stats_2x3 <- function(data, outcome, group, pcutoff, padj.method,
                           test = c("t.test", "wilcoxon"),
                           exact = FALSE) {
  test <- match.arg(test)
  res <- lapply(levels(group), function(i) {
    subdat <- data[group == i, ]
    suboutcome <- as.numeric(outcome[group == i])
    out <- col_t_welch(subdat[suboutcome==1, ], subdat[suboutcome==2, ])
    if (test == "wilcoxon") {
      out$pvalue <- col_wilcoxon_twosample(subdat[suboutcome==1, ],
                                           subdat[suboutcome==2, ],
                                           exact = exact)$pvalue
    }
    out$padj <- qval(out$pvalue, padj.method)
    out
  })
  names(res) <- levels(group)
  res
}
