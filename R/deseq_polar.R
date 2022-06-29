#' Convert DESeq2 objects to a volcano3d object
#'
#' This function takes 2 `DESeqDataSet` objects and converts the results to a
#' 'volc3d' object.
#'
#' @param object An object of class 'DESeqDataSet' with the full design formula.
#'   The function `DESeq` needs to have been run.
#' @param objectLRT An object of class 'DESeqDataSet' with the reduced design
#'   formula. The function `DESeq` needs to have been run on this object with
#'   argument `test="LRT"`.
#' @param contrast Character value specifying column within the metadata stored
#'   in the DESeq2 dataset objects is the outcome variable. This column must 
#'   contain a factor with 3 levels.
#' @param data Optional matrix containing gene expression data. If not supplied,
#'   the function will pull the expression data from within the DESeq2 object
#'   using the DESeq2 function `assay()`. NOTE: for consistency with gene
#'   expression datasets, genes are in rows.
#' @param pcutoff Cut-off for p-value significance
#' @param padj.method Can be any method available in `p.adjust` or `"qvalue"`.
#'   The option "none" is a pass-through.
#' @param ... Optional arguments passed to [polar_coords()]
#' @export

deseq_polar <- function(object, objectLRT, contrast,
                        data = NULL,
                        pcutoff = 0.05,
                        padj.method = "BH", ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Can't find package DESeq2. Try:
           BiocManager::install('DESeq2')",
         call. = FALSE)
  }
  if (!inherits(object, "DESeqDataSet")) stop("Not a DESeqDataSet object")
  if (!inherits(objectLRT, "DESeqDataSet")) stop("Not a DESeqDataSet object")
  LRT <- DESeq2::results(objectLRT)
  LRT <- as.data.frame(LRT[, c('pvalue', 'padj')])
  if (padj.method == "qvalue") {
    LRT <- deseq_qvalue(LRT)
    index <- LRT$qvalue < pcutoff & !is.na(LRT$qvalue)
    ptype <- "qvalue"
  } else {
    index <- LRT$padj < pcutoff & !is.na(LRT$padj)
    ptype <- "padj"
  }
  groups <- levels(object@colData[, contrast])
  contrastlist <- list(
    groups[1:2],
    groups[c(3, 1)],
    groups[2:3]
  )
  pairres <- lapply(contrastlist, function(i) {
    res <- DESeq2::results(object, contrast = c(contrast, i))
    as.data.frame(res[, c('pvalue', 'padj')])
  })
  pvals <- cbind(LRT[, "pvalue"], pairres[[1]][, "pvalue"], 
                 pairres[[2]][, "pvalue"], pairres[[3]][, "pvalue"])
  pairadj <- lapply(pairres, function(res) {
    out <- rep_len(NA, nrow(LRT))
    out[index] <- if (padj.method == "qvalue") {
      deseq_qvalue(res[index, ])$qvalue
    } else p.adjust(res[index, 'pvalue'], method = padj.method)
    out
  })
  padj <- cbind(LRT[, ptype], pairadj[[1]], pairadj[[2]], pairadj[[3]])
  dimnames(pvals) <- dimnames(padj) <- list(rownames(LRT), c("LRT", "AvB", "AvC", "BvC"))
  if (is.null(data)) {
    vstdata <- DESeq2::vst(object)
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Can't find package SummarizedExperiment. Try:
           BiocManager::install('SummarizedExperiment')",
           call. = FALSE)
    }
    data <- SummarizedExperiment::assay(vstdata)
  }
  polar_coords(object@colData[, contrast], t(data), pvals, padj, pcutoff, ...)
}


deseq_qvalue <- function(df) {
  q <- qval(df$pvalue[!is.na(df$padj)])
  df$qvalue <- 1  # NA converted to 1
  df$qvalue[!is.na(df$padj)] <- q
  df
}


deseq_LRT <- function(object, formula, ...) {
  
}
