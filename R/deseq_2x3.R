
#' 2 x 3 factor DESeq2 analysis
#' 
#' Experimental function for performing 2x3 factor DESeq2 analyses. Output can
#' be passed to [deseq_2x3_polar()] and subsequently plotted. Example usage
#' would include comparing gene expression against a binary outcome e.g.
#' response vs non-response, across 3 drugs: the design would be `~ response`
#' and `group` would refer to the medication column in the metadata.
#' 
#' @param object An object of class 'DESeqDataSet' containing full dataset
#' @param design Design formula. The main contrast is taken from the last term
#'   of the formula and must be a binary factor.
#' @param group Character value for the column with the 3-way grouping factor
#'   within the sample information data `colData`
#' @param ... Optional arguments passed to `DESeq()`.
#' @return Returns a list of 3 DESeq2 results objects which can be passed onto
#'   [deseq_2x3_polar()]
#' @examples
#' \donttest{
#' # Basic DESeq2 set up
#' library(DESeq2)
#' 
#' counts <- matrix(rnbinom(n=3000, mu=100, size=1/0.5), ncol=30)
#' rownames(counts) <- paste0("gene", 1:100)
#' cond <- rep(factor(rep(1:3, each=5), labels = c('A', 'B', 'C')), 2)
#' resp <- factor(rep(1:2, each=15), labels = c('non.responder', 'responder'))
#' metadata <- data.frame(drug = cond, response = resp)
#' 
#' # Full dataset object construction
#' dds <- DESeqDataSetFromMatrix(counts, metadata, ~response)
#' 
#' # Perform 3x DESeq2 analyses comparing binary response for each drug
#' res <- deseq_2x3(dds, ~response, "drug")
#' 
#' # Generate polar object
#' obj <- deseq_2x3_polar(res)
#' 
#' # 2d plot
#' radial_plotly(obj)
#' 
#' # 3d plot
#' volcano3D(obj)
#' }
#' 
#' @export

deseq_2x3 <- function(object, design, group, ...) {
  # Check packages and input data
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Can't find package DESeq2. Try:
           BiocManager::install('DESeq2')", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Can't find package SummarizedExperiment. Try:
           BiocManager::install('SummarizedExperiment')",
         call. = FALSE)
  }
  if (!inherits(object, "DESeqDataSet")) stop("Not a DESeqDataSet object")
  counts <- SummarizedExperiment::assay(object)
  colDat <- SummarizedExperiment::colData(object)
  termsDE <- attr(terms(design), "term.labels")
  contrast <- termsDE[length(termsDE)]
  if (nlevels(colDat[, contrast]) != 2) stop(contrast, " is not binary in `design`")
  if (!group %in% colnames(colDat)) {
    stop(group, " is not a column in sample information in `colData`")}
  groups <- colDat[, group]
  if (nlevels(groups) != 3) stop(group, " does not have 3 levels")
  res <- lapply(levels(groups), function(i) {
    message(group, " = ", i)
    dds <- DESeq2::DESeqDataSetFromMatrix(counts[, groups == i],
                                          colDat[groups == i, ], design)
    dds <- DESeq2::DESeq(dds, ...)
    DESeq2::results(dds)
  })
  setNames(res, levels(groups))
}
