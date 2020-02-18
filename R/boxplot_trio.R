#' Boxplot to compare three groups
#'
#' Plots the expression of a specific row in expression comparing the three 
#' groups in a boxplot. 
#' @param dep A dep object with the pvalues between groups of interest. Created
#' by \code{\link{create_dep}}.
#' @param value The row name or number in dep@@expression to be analysed
#' @param box_colours The fill colours for each box 
#' (default = c('blue', 'red', 'green3') ).
#' @param test The statistical test used to compare expression. See method in
#' \code{\link[ggpubr]{stat_compare_means}}
#' @param levels_order A character vector stating the contrast groups to be 
#' plotted in order
#' @param my_comparisons A list of contrasts to pass to comparisons in 
#' \code{\link[ggpubr]{stat_compare_means}}. If NULL (default) all contrast 
#' pvalues are calculated and plotted. 
#' @param text_size The font size of text (default = 10)
#' @param stat_size The font size of statistical parameter (default = 3).
#' @param ... Other parameters for \code{\link[ggpubr]{stat_compare_means}}
#' @return Returns a ggplot boxplot featuring the expression of one row in 
#' expression across the three groups in comparison with annotated pvalues. 
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @importFrom ggplot2 theme ggplot labs geom_path geom_path geom_text annotate 
#' geom_point scale_color_manual aes geom_jitter
#' @importFrom utils combn
#' @keywords pvalue, polar, plot, boxplot
#' @references 
#' Lewis, Myles J., et al. (2019). 
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and 
#' treatment response phenotypes.} 
#' \emph{Cell reports}, \strong{28}:9 
#' @export
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
#' 
#' boxplot_trio(syn_p_obj,
#'           value = "SLAMF6",
#'           test = "wilcox.test",
#'           levels_order = c("Lymphoid", "Myeloid", "Fibroid"),
#'           box_colours = c("blue", "red", "green3"))
#' 
#' boxplot_trio(syn_p_obj,
#'           value = "ITM2C",
#'           test = "wilcox.test",
#'           levels_order = c("Lymphoid", "Myeloid", "Fibroid"),
#'           my_comparisons = list(c("Lymphoid", "Myeloid"),
#'                                 c("Myeloid", "Fibroid")),
#'           box_colours = c("blue", "red", "green3"))


boxplot_trio <- function(dep, 
                         value, 
                         box_colours = c('green3', 'blue', 'red'), 
                         test = "t.test", 
                         levels_order = NULL, 
                         my_comparisons = NULL,
                         text_size = 10,
                         stat_size = 3, 
                         ...){
    
sampledata <- dep@sampledata
    expression <- dep@expression
    
    if(is.null(levels_order)) levels_order <- levels(sampledata[, dep@contrast])
    
    if(! class(levels_order) %in% c("character")) {
        stop("sampledata must be a data frame")
    }
    if(! all(levels_order %in% levels(sampledata[, dep@contrast]))){
        stop('levels_order must be a character vector defining the order of 
        levels in sampledata[, contrast]')
    }
    if(! class(sampledata) %in% c("data.frame")) {
        stop("sampledata must be a data frame")
    }
    if(! class(expression) %in% c("data.frame", "matrix")) {
        stop("expression must be a data frame or matrix")
    }
    if(! class(value) %in% c("character", "numeric")) {
        stop("value must be a character")
    }
    if(length(value) > 1) stop("value must be of length 1")
    if(! value %in% rownames(expression)) {
        stop("value/gene is not in rownames(expression)")
    }
    if(! identical(colnames(expression), as.character(sampledata$ID))) {
        stop("expression and sampledata misalligned")
    }
    
    sampledata$comp <- sampledata[, dep@contrast]
    sampledata <- sampledata[sampledata$comp %in% levels_order, ]
    sampledata$comp <- factor(sampledata$comp, levels = levels_order)
    
    expression <- expression[, match(as.character(sampledata$ID), 
                                     colnames(expression))]
    
    if(class(value) ==  "character") {
        index <- which(rownames(expression) ==  value)
    }
    
    if(is.null(my_comparisons)) {
        comps <- levels(sampledata$comp)
        my_comparisons <- lapply(seq_len(ncol(combn(comps, 2))), function(i) {
            as.character(combn(comps, 2)[, i])
        })
    }
    df <- data.frame("group" = sampledata$comp, 
                     "row" = as.numeric(as.character(expression[value, ])))
    df <- df[! is.na(df$row), ]
    
    ggboxplot(data = df, 
              x = "group", 
              y = "row", 
              xlab = "", 
              ylab = rownames(expression)[index],
              fill = "group", 
              palette = box_colours, 
              outlier.shape = NA, 
              alpha = 0.5) +
        geom_jitter(height = 0, width = 0.30) +
        stat_compare_means(comparisons = my_comparisons, method = test, 
                           size=stat_size, ...) +
        theme(legend.position = "none", 
              text = element_text(size = text_size))
}
