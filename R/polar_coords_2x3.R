

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
  
  res <- calc_pvals_2x3(data, outcome, group, pcutoff, padj.method, ...)
  
  df1 <- getCols(res, rn, 'statistic')
  df2 <- getCols(res, rn, 'mean.diff')
  pvals <- getCols(res, rn, 'pvalue')
  padj <- getCols(res, rn, 'padj')
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
               df = list(scaled = df1, unscaled = df2, type = "2x3_polar"),
               outcome = outcome,
               data = data.frame(), pvals = pvals, padj = padj,
               pcutoff = pcutoff, scheme = scheme,
               labs = levels(ptab$lab))
}



#' @importFrom matrixTests col_t_welch col_wilcoxon_twosample
#' @export
#' 
calc_pvals_2x3 <- function(data, outcome, group, pcutoff, padj.method,
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
