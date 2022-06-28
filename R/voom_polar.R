

#' Convert RNA-Seq count data to a volcano3d object using 'limma voom'
#'
#' This function takes a design formula, metadata and raw count data and uses
#' 'limma voom' to analyse the data. The results are converted to a 'volc3d'
#' object ready for plotting a 3d volcano plot or polar plot.
#'
#' @param formula Design formula which must be of the form `~ 0 + outcome + ...`.
#'   The 3-way outcome variable must be the first variable after the '0', and
#'   this variable must be a factor with exactly 3 levels.
#' @param metadata Matrix or dataframe containing metadata as referenced by
#'   `formula`
#' @param counts Matrix containing raw gene expression count data
#' @param pcutoff Cut-off for p-value significance
#' @param padj.method Can be any method available in `p.adjust` or `"qvalue"`.
#'   The option "none" is a pass-through.
#' @param ... Optional arguments passed to [polar_coords()]
#' @export

voom_polar <- function(formula, metadata, counts,
                       pcutoff = 0.05,
                       padj.method = "BH", ...) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Can't find package edgeR. Try:
           BiocManager::install('edgeR')", call. = FALSE)
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Can't find package limma. Try:
           BiocManager::install('limma')", call. = FALSE)
  }
  modterms <- attr(terms(formula), "term.labels")
  outcome_col <- modterms[1]
  if (nlevels(metadata[, outcome_col]) != 3) stop("Outcome does not have 3 levels")
  
  design <<- model.matrix(formula, data = metadata)  # needed for limma
  contrast.matrix <- limma::makeContrasts(
    paste0(colnames(design)[1] , "-", colnames(design)[2]),
    paste0(colnames(design)[1] , "-", colnames(design)[3]),
    paste0(colnames(design)[2] , "-", colnames(design)[3]),
    levels = design)
  dge <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  # voom
  v <- limma::voom(dge, design, plot = FALSE)
  fit1 <- limma::lmFit(v, design)
  fit <- limma::contrasts.fit(fit1, contrast.matrix) 
  fit <- limma::eBayes(fit)
  contrasts <- colnames(coefficients(fit))
  Pvals_limma_DE <- lapply(contrasts, function(x){
    id <- which(colnames(coefficients(fit)) == x)
    out <- limma::topTable(fit, coef= id, number=Inf, sort.by="none")
    out[,c("P.Value", "logFC")]
  })
  out <- limma::topTable(fit, coef=1:3, number=Inf, sort.by="none")
  Pvals_overall <- out$P.Value
  pvals <- cbind(Pvals_overall,
                 Pvals_limma_DE[[1]]$P.Value, 
                 Pvals_limma_DE[[2]]$P.Value, 
                 Pvals_limma_DE[[3]]$P.Value)
  LRTpadj <- qval(Pvals_overall, method = padj.method)
  ind <- LRTpadj < pcutoff
  pairadj <- apply(pvals[, 2:4], 2, function(res) {
    out <- rep_len(NA, length(LRTpadj))
    out[ind] <- qval(res[ind], method = padj.method)
    out
  })
  padj <- cbind(LRTpadj, pairadj)
  
  polar_coords(metadata[, outcome_col], t(log2(counts[keep, ] + 1)), pvals, padj, pcutoff, ...)
}

