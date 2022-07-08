setClassUnion("df_or_matrix", c("data.frame", "matrix"))

setClass("volc3d", slots = list(df = "list",
                                outcome = "factor",
                                data = "df_or_matrix",
                                pvals = "matrix",
                                padj = "df_or_matrix",
                                pcutoff = "numeric",
                                scheme = "character",
                                labs = "character"))


#' Coordinates for Three Way Polar Plot
#'
#' This function creates a 'volc3d' object of S3 class for downstream plots
#' containing the p-values from a three-way group comparison, expression data
#' sample data and polar coordinates. For RNA-Seq count data, two functions
#' \code{\link{deseq_polar}} or \code{\link{voom_polar}} can be used instead.
#'
#' @param outcome Outcome vector with 3 groups, ideally as a factor. If it is
#'   not a factor, this will be coerced to a factor. This must have exactly 3
#'   levels. NOTE: if `pvals` is given, the order of the levels in `outcome`
#'   must correspond to the order of columns in `pvals`.
#' @param data Dataframe or matrix with variables in columns
#' @param pvals Matrix or dataframe with p-values. The first column represents a
#'   test across all 3 categories such as one-way ANOVA or likelihood ratio
#'   test. Columns 2-4 represent pairwise tests comparing groups A vs B, A vs C
#'   and B vs C, where A, B, C represent levels 1, 2, 3 in `outcome`. Columns
#'   2-4 must be provided in the correct order. If `pvals` is not given, it is
#'   calculated using the function \code{\link{calc_pvals}}.
#' @param padj Matrix or dataframe with adjusted p-values. If not supplied,
#'   defaults to use nominal p-values from `pvals`.
#' @param pcutoff Cut-off for p-value significance
#' @param scheme Vector of colours starting with non-significant variables
#' @param labs Optional character vector for labelling groups. Default `NULL`
#'   leads to abbreviated labels based on levels in `outcome` using
#'   [abbreviate()]. A vector of length 3 with custom abbreviated names for the
#'   outcome levels can be supplied. Otherwise a vector length 7 is expected, of
#'   the form "ns", "B+", "B+C+", "C+", "A+C+", "A+", "A+B+", where "ns" means
#'   non-significant and A, B, C refer to levels 1, 2, 3 in `outcome`, and must
#'   be in the correct order.
#' @param ... Optional arguments passed to \code{\link{calc_pvals}}
#' 
#' @return Returns an S4 'volc3d' object containing:
#' \itemize{
#'   \item{'df'} A list of 2 dataframes. Each dataframe contains both x,y,z
#'   coordinates as well as polar coordinates r, angle. The first dataframe has
#'   coordinates on scaled data. The 2nd dataframe has unscaled data (e.g. log2
#'   fold change for gene expression).
#'   \item{'outcome'} The three-group contrast factor used for comparisons
#'   \item{'data'} Dataframe or matrix containing the expression data
#'   \item{'pvals'} A dataframe containing p-values. First column is the 3-way
#'   comparison (LRT or ANOVA). Columns 2-4 are pairwise comparisons between
#'   groups A vs B, A vs C and B vs C, where A, B, C are the 3 levels in the
#'   outcome factor.
#'   \item{'padj'} A dataframe containing p-values adjusted for multiple testing
#' }
#' 
#' @seealso \code{\link{deseq_polar}}, \code{\link{voom_polar}}
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
#'                           data = t(syn_example_rld))
#' 
#' @importFrom methods is
#' @export
#'
polar_coords <- function(outcome, data,
                       pvals = NULL, padj = pvals, pcutoff = 0.05,
                       scheme = c('grey60', 'red', 'gold2', 'green3', 
                                  'cyan', 'blue', 'purple'),
                       labs = NULL, ...) {
  if (length(outcome) != nrow(data)) {
    stop("Number of rows in `data` differs from `outcome`")}
  if (any(is.na(outcome))) {
    ok <- !is.na(outcome)
    data <- data[ok,]
    outcome <- outcome[ok]
    message("Removing NA from `outcome`")
  }
  outcome <- as.factor(outcome)
  outcome <- droplevels(outcome)
  if (nlevels(outcome) != 3) stop("`outcome` must have 3 levels")
  data <- as.matrix(data)
  data_sc <- scale(data)
  df1 <- vapply(levels(outcome), function(i) colMeans(data_sc[outcome == i, ]),
                numeric(ncol(data)))
  df2 <- vapply(levels(outcome), function(i) colMeans(data[outcome == i, ]),
                numeric(ncol(data)))
  df1 <- polar_xy(df1)
  df2 <- polar_xy(df2)
  if (is.null(pvals)) {
    pv <- calc_pvals(outcome, data, pcutoff, ...)
    pvals <- pv$pvals
    padj <- pv$padj
  }
  ptab <- polar_p(outcome, df1, pvals, padj, pcutoff, scheme, labs)
  df1 <- cbind(df1, ptab)
  df2 <- cbind(df2, ptab)
  methods::new("volc3d",
               df = list(scaled = df1, unscaled = df2),
               outcome = outcome, data = data, pvals = pvals, padj = padj,
               pcutoff = pcutoff, scheme = scheme,
               labs = levels(ptab$lab))
}


polar_xy <- function(df, angle_offset = 0) {
  y <- sinpi(1/3) * (df[,2] - df[,3])
  x <- df[,1] - (cospi(1/3) * (df[,3] + df[,2]))
  r <- sqrt(x^2 + y^2)
  angle <- atan2(y, x)/(2*pi)
  angle <- ((angle + angle_offset) %% 1) * 360
  cbind(df, x, y, r, angle)
}


#' Calculate one-way test and pairwise tests
#' 
#' Internal function for calculating 3-class group test (either one-way ANOVA or
#' Kruskal-Wallis test) and pairwise tests (either t-test or Wilcoxon test) on
#' multi-column data against an outcome parameter with 3 levels.
#' 
#' @param outcome Outcome vector with 3 groups, ideally as a factor. If it is
#'   not a factor, this will be coerced to a factor. This must have exactly 3
#'   levels.
#' @param data Dataframe or matrix with variables in columns
#' @param pcutoff Cut-off for p-value significance
#' @param padj.method Can be any method available in `p.adjust` or `"qvalue"`.
#'   The option "none" is a pass-through.
#' @param group_test Specifies statistical test for 3-class group comparison.
#'   "anova" means one-way ANOVA, "kruskal.test" means Kruskal-Wallis test.
#' @param pairwise_test Specifies statistical test for pairwise comparisons
#' @param exact Logical which is only used with `pairwise_test = "wilcoxon"`
#' @param filter_pairwise Logical. If `TRUE` (the default) p-value adjustment on
#'   pairwise statistical tests is only conducted on attributes which reached
#'   the threshold for significance after p-value adjustment on the group
#'   statistical test.
#' @importFrom Rfast ftests ttests kruskaltests
#' @importFrom matrixTests row_wilcoxon_twosample
#' @export
#'
calc_pvals <- function(outcome, data,
                       pcutoff = 0.05,
                       padj.method = "BH",
                       group_test = c("anova", "kruskal.test"),
                       pairwise_test = c("t.test", "wilcoxon"),
                       exact = FALSE,
                       filter_pairwise = TRUE) {
  group_test <- match.arg(group_test)
  pairwise_test <- match.arg(pairwise_test)
  outcome <- as.factor(outcome)
  if (nlevels(outcome) != 3) stop("`outcome` must have 3 levels")
  data <- as.matrix(data)
  res <- switch(group_test,
                "anova" = Rfast::ftests(data, outcome),
                "kruskal.test" = Rfast::kruskaltests(data, outcome))
  rownames(res) <- colnames(data)
  onewayp <- res[, 2]
  indx <- lapply(levels(outcome), function(i) outcome == i)
  if (pairwise_test == "wilcoxon") {
    res1 <- suppressWarnings(
      matrixTests::row_wilcoxon_twosample(t(data[indx[[1]], ]), t(data[indx[[2]], ]),
                                          exact = exact))
    res2 <- suppressWarnings(
      matrixTests::row_wilcoxon_twosample(t(data[indx[[1]], ]), t(data[indx[[3]], ]),
                                          exact = exact))
    res3 <- suppressWarnings(
      matrixTests::row_wilcoxon_twosample(t(data[indx[[2]], ]), t(data[indx[[3]], ]),
                                          exact = exact))
  } else {
    res1 <- Rfast::ttests(data[indx[[1]], ], data[indx[[2]], ])
    res2 <- Rfast::ttests(data[indx[[1]], ], data[indx[[3]], ])
    res3 <- Rfast::ttests(data[indx[[2]], ], data[indx[[3]], ])
  }
  p1 <- res1[, "pvalue"]
  p2 <- res2[, "pvalue"]
  p3 <- res3[, "pvalue"]
  pvals <- cbind(onewayp, p1, p2, p3)
  if (padj.method == "none") {
    padj <- pvals
  } else {
    onewaypadj <- qval(onewayp, method = padj.method)
    index <- if (filter_pairwise) {
      onewaypadj < pcutoff & !is.na(onewaypadj)
    } else !is.na(onewaypadj)
    padj <- data.frame(onewaypadj, p1 = NA, p2 = NA, p3 = NA)
    padj$p1[index] <- qval(p1[index], method = padj.method)
    padj$p2[index] <- qval(p2[index], method = padj.method)
    padj$p3[index] <- qval(p3[index], method = padj.method)
  }
  padj <- as.matrix(padj)
  
  list(pvals = pvals, padj = padj)
}


#' @importFrom stats p.adjust p.adjust.methods
#'
qval <- function(p, method = "qvalue") {
  if (method %in% p.adjust.methods) return(p.adjust(p, method = method))
  if (!requireNamespace("qvalue", quietly = TRUE)) {
    stop("Can't find package qvalue. Try:
           BiocManager::install('qvalue')",
         call. = FALSE)
  }
  q <- try(qvalue::qvalue(p)$qvalues, silent = TRUE)
  if (inherits(q, 'try-error')) q <- p.adjust(p, method = "BH")
  q
}


#' @importFrom Rfast rowMins
#'
polar_p <- function(outcome, df1, pvals, padj = pvals, pcutoff = 0.05,
                    scheme = c('grey60', 'red', 'gold2', 'green3', 
                               'cyan', 'blue', 'purple'),
                    labs = NULL) {
  pvalue <- pvals[,1]
  z <- -log10(pvals[,1])
  paircut <- padj[, 2:4] < pcutoff
  paircut <- paircut *1  # convert matrix to numeric
  mincol <- Rfast::rowMins(as.matrix(df1[, 1:3]))
  mincol2 <- c("A", "B", "C")[mincol]
  pairmerge <- paste0(mincol2, paircut[,1], paircut[,2], paircut[,3])
  pgroup <- rep_len(1, nrow(df1))
  # sequence AB, AC, BC
  pgroup[grep("A10.|C.01", pairmerge)] <- 2  # red
  pgroup[grep("A01.|B0.1", pairmerge)] <- 4  # green
  pgroup[grep("B1.0|C.10", pairmerge)] <- 6  # blue
  pgroup[grep("A11.", pairmerge)] <- 3  # yellow
  pgroup[grep("B1.1", pairmerge)] <- 5  # cyan
  pgroup[grep("C.11", pairmerge)] <- 7  # purple
  pgroup[pvals[,1] > pcutoff] <- 1  # ns for all p_lrt < cutoff
  col <- scheme[pgroup]
  if (is.null(labs) | length(labs) == 3) {
    abbrev <- if (length(labs) == 3) labs else abbreviate(levels(outcome), 1)
    labs <- c("ns",
              paste0(abbrev[2], "+"),
              paste0(abbrev[2], "+", abbrev[3], "+"),
              paste0(abbrev[3], "+"),
              paste0(abbrev[1], "+", abbrev[3], "+"),
              paste0(abbrev[1], "+"),
              paste0(abbrev[1], "+", abbrev[2], "+"))
  }
  lab <- factor(pgroup, levels = 1:7, labels = labs)
  data.frame(z, pvalue, col, lab)
}
