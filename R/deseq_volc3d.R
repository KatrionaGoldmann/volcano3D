#' Convert DESeq2 objects to volcano3d
#'
#'
#'
#' @importFrom DESeq2 results vst
#' @importFrom SummarizedExperiment assay
#' @export

DESeqToVolc <- function(object, objectLRT, contrast,
                        data = NULL,
                        pcutoff = 0.05,
                        padj.method = "BH", ...) {
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
    res <- results(object, contrast = c(contrast, i))
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
    data <- SummarizedExperiment::assay(vstdata)
  }
  polarCoord(object@colData[, contrast], t(data), pvals, padj, pcutoff, ...)
}


deseq_qvalue <- function(df) {
  q <- qval(df$pvalue[!is.na(df$padj)])
  df$qvalue <- 1  # NA converted to 1
  df$qvalue[!is.na(df$padj)] <- q
  df
}


#' @export
quick_volcano3d <- function(obj, type = 1) {
  if (!inherits(obj, "volc3d")) stop("Not a 'volc3d' class object")
  plot_ly(obj[[type]], x = ~x, y = ~y, z = ~z, color = ~lab, colors = obj$scheme,
          hoverinfo='text',
          text = ~paste0(rownames(obj[[type]]), "<br>theta = ", as.integer(angle),
                         "<br>r = ", formatC(r, digits = 3),
                         "<br>P = ", format(pvalue, digits = 3, scientific = 3)),
          marker = list(size = 3),
          type = "scatter3d", mode = "markers")
}
