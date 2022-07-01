

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
#' @param filter_pairwise Logical whether adjusted p-value pairwise statistical
#'   tests are only conducted on genes which reach significant adjusted p-value
#'   cut-off on the group likelihood ratio test
#' @param ... Optional arguments passed to [polar_coords()]
#' 
#' @examples
#' if (requireNamespace("limma", quietly = TRUE) & 
#'     requireNamespace("edgeR", quietly = TRUE)) {
#'   library(limma)
#'   library(edgeR)
#'   
#'   counts <- matrix(rnbinom(n=1500, mu=100, size=1/0.5), ncol=15)
#'   cond <- factor(rep(1:3, each=5), labels = c('A', 'B', 'C'))
#'   cond <- data.frame(cond)
#'   
#'   polar <- voom_polar(~0 + cond, cond, counts)
#'   
#'   volcano3D(polar)
#'   radial_ggplot(polar)
#' }
#' 
#' @importFrom stats coefficients model.matrix terms update
#' @export

voom_polar <- function(formula, metadata, counts,
                       pcutoff = 0.05,
                       padj.method = "BH",
                       filter_pairwise = TRUE, ...) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Can't find package edgeR. Try:
           BiocManager::install('edgeR')", call. = FALSE)
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Can't find package limma. Try:
           BiocManager::install('limma')", call. = FALSE)
  }
  # force formula to have no intercept
  if (attr(terms(formula), "intercept") != 0) formula <- update(formula, ~ . -1)
  modterms <- attr(terms(formula), "term.labels")
  outcome_col <- modterms[1]
  if (nlevels(metadata[, outcome_col]) != 3) stop("Outcome does not have 3 levels")
  
  vdesign <- model.matrix(formula, data = metadata)
  contrast_set <- list(paste0(colnames(vdesign)[1] , "-", colnames(vdesign)[2]),
                       paste0(colnames(vdesign)[1] , "-", colnames(vdesign)[3]),
                       paste0(colnames(vdesign)[2] , "-", colnames(vdesign)[3]))
  if(!exists("makeContrasts")) stop("limma package not loaded")
  contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),
                                    contrast_set, levels=list(vdesign))))
  
  dge <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, vdesign)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  # voom
  v <- limma::voom(dge, vdesign, plot = FALSE)
  fit1 <- limma::lmFit(v, vdesign)
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
  ind <- if (filter_pairwise) LRTpadj < pcutoff else rep_len(TRUE, length(LRTpadj))
  pairadj <- apply(pvals[, 2:4], 2, function(res) {
    out <- rep_len(NA, length(LRTpadj))
    out[ind] <- qval(res[ind], method = padj.method)
    out
  })
  padj <- cbind(LRTpadj, pairadj)
  
  polar_coords(metadata[, outcome_col], t(log2(counts[keep, ] + 1)), pvals, padj, pcutoff, ...)
}

