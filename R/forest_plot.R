
#' Forest plot individual gene from 2x3 factor analysis
#' 
#' @param object A 'volc3d' class object
#' @param gene Gene to plot
#' @param scheme Vector of 3 colours for plotting
#' @param ... Optional arguments passed to [plot()]
#' @importFrom graphics abline arrows axis par text
#' @export

forest_plot <- function(object, gene,
                        scheme = c('red', 'green3', 'blue'),
                        ...) {
  df <- object@df[[2]]
  x <- unlist(df[gene, 1:3])
  xCI <- unlist(df[gene, 4:6]) * 1.96
  xrange <- range(c(x-xCI, x+xCI, 0), na.rm=TRUE)
  op <- par(mar=c(5, 7, 2, 4))
  on.exit(par(op))
  plot(x, y = 3:1,
       pch=19, col=scheme, 
       xlim = xrange, yaxt="n", xlab=paste(gene, "log2 FC"), ylab = "", las=1, bty="n",
       panel.first = {
         arrows(x-xCI, 3:1, x+xCI, code=3, angle=90, length= 0.05,
                col=scheme)
       }, ...)
  abline(v=0, lty=2)
  axis(2, 3:1, colnames(df)[1:3], tick=F, las=1)
  stars <- stars.pval(object@padj[gene,])
  stars[stars== "."] <- ""
  text(xrange[2]+ diff(xrange)*0.06, 3:1, stars, xpd = NA)
}


# code modified from orphaned package gtools

#' @importFrom stats symnum
stars.pval <- function(p.value) {
  unclass(
    symnum(p.value,
           corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 1),
           symbols = c("***", "**", "*", " ")
    )
  )
}
