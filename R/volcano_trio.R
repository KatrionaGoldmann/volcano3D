#' Creates a single volcano plot
#'
#' This function creates a volcano plot for one comparison group
#' @param pvalues_df The pvalues data frame. This must contain pvalue, padj, 
#' logFC columns as well as a label column. 
#' @param comparison The comparison to use
#' @importFrom ggplot2 ggplot labs geom_path geom_path geom_text annotate 
#' geom_point scale_colour_manual theme aes theme_classic element_text 
#' geom_hline geom_vline unit layer_scales lims element_rect aes_string
#' @param p_cutoff The cut-off for pvalue significance (default = 0.05). 
#' @param fc_cutoff The cut-off for fold change significance (default = 1). 
#' @param p_col_suffix The suffix word to define columns containing p values
#' (default = 'pvalues').
#' @param padj_col_suffix The suffix word to define columns containing adjusted
#' p values (default = 'padj'). If NULL these will be calculated using
#' \code{padjust_method}.
#' @param fc_col_suffix The optional suffix word to define columns containing 
#' log fold change values (default = 'logFC').
#' @param label_column Optional column name in pvalues for markers to be 
#' labelled with at plotting stage. If NULL the rownames of pvalues are used. 
#' @param label_rows Row numbers or names for values to be annotated/labelled
#' (default = NULL).
#' @param label_size The font size of labels (default = 3)
#' @param text_size The font size of text (default = 10)
#' @param marker_size The size of markers (default = 0.7)
#' @param shared_legend_size The size for the legend (default = 1). 
#' @param sig_names A character vector of labels to be used for: 
#' non-significant; adjusted p < p_cutoff; |Fold Change| > fc_cutoff; and 
#' finally adjusted p < p_cutoff. 
#' If NULL c('Not Significant', \code{paste('Padj <', `p_cutoff`)},
#' \code{paste('|FC| >', `fc_cutoff`)},
#' \code{paste('Padj <', `p_cutoff`, 'and |FC| >', `fc_cutoff`)}) is used.
#' @param colours A character vector of colours to be used for non-significant;
#' adjusted p < p_cutoff;
#' |Fold Change| >, fc_cutoff; and adjusted p < p_cutoff. default = p_cutoff &
#' |Fold Change| > fc_cutoff markers respectively
#' (default = c('grey60', 'salmon', 'steelblue', 'limegreen')).
#' @param fc_line Logical whether to add vertical dashed line at fc_cutoff
#' (default = TRUE).
#' @param p_line Logical whether to add horizontal dashed line at p_cutoff
#' (default = TRUE).
#' @param line_colours A character vector stating the colour of lines to be
#' used for fc_line and p_line respectively (default = c('black', 'black')).
#' @importFrom ggplot2 ggplot labs geom_path geom_path geom_text annotate 
#' geom_point scale_colour_manual theme aes theme_classic element_text 
#' geom_hline geom_vline unit element_rect aes_string
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats p.adjust setNames
#' @return Returns a single ggplot volcano plots. 
#' @concept volcanoplot 
#' @keywords hplot
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @export
#' @examples 
#' data("example_data")
#' volcano_plot(syn_example_p, 
#'              "Fibroid_Lymphoid", 
#'              label_col = "Gene", 
#'              label_rows=c("SLAMF6"), 
#'              fc_col_suffix="log2FoldChange") 

volcano_plot <- function(pvalues_df, 
                         comparison, 
                         p_cutoff = 0.05,
                         fc_cutoff = 1,
                         label_rows = NULL,
                         p_col_suffix = "pvalue",
                         padj_col_suffix = "padj",
                         fc_col_suffix = "logFC",
                         label_col = "label",
                         label_size = 3,
                         text_size = 10,
                         marker_size = 0.7,
                         shared_legend_size = 1, 
                         sig_names = NULL,
                         colours = c("grey60", "salmon", 
                                     "steelblue", "limegreen"),
                         fc_line = TRUE,
                         p_line = TRUE,
                         line_colours = c("black", "black")){
  
  toptable <- pvalues_df[, c(
    colnames(pvalues_df)[grepl(comparison, colnames(pvalues_df))], label_col)]
  colnames(toptable)[colnames(toptable) == label_col] <- "label"
  
  if(is.null(label_col)){
    toptable$label <- rownames(toptable)
  }
  
  if(is.null(sig_names)){
    sig_names <- c("Not Significant",
                   paste("Padj <", p_cutoff),
                   paste("|log(FC)| >", fc_cutoff),
                   paste("Padj <", p_cutoff, "and |log(FC)| >", fc_cutoff))
  }
  
  colnames(toptable)[unlist(
    lapply(c(p_col_suffix, padj_col_suffix, fc_col_suffix), function(x) {
      which(grepl(x, colnames(toptable)))
    }))] <- c("pvalue", "padj", "logFC")
  
  toptable$cols <- sig_names[1]
  toptable$cols[toptable$padj <=  p_cutoff] <- sig_names[2]
  toptable$cols[abs(toptable$logFC) > fc_cutoff] <- sig_names[3]
  toptable$cols[toptable$padj <= p_cutoff &
                  abs(toptable$logFC) > fc_cutoff] <- sig_names[4]
  toptable$cols <- factor(toptable$cols)
  mapping <- setNames(sig_names, colours)
  
  toptable <- toptable[! is.na(toptable$logFC), ]
  toptable <- toptable[! is.na(-log10(toptable$pvalue)), ]
  toptable <- toptable[! is.na(-log10(toptable$padj)), ]
  toptable$lp <- -log10(toptable$pvalue)
  
  # Create the volcano plot ggplot
  p <- ggplot(toptable, 
              aes_string(x = "logFC", 
                         y = "lp"),
              colour = "cols") +
    geom_point(data = toptable, size=marker_size, aes_string(color="cols")) +
    scale_color_manual(values = names(mapping)[match(levels(toptable$cols),
                                                     mapping)]) +
    labs(y = expression(-log["10"]*p), 
         x = expression(log["2"]*FC),
         title = gsub("_", " vs ", comparison), 
         color = "Significance") +
    theme_classic() +
    theme(text = element_text(size = text_size), 
          legend.background = element_rect(fill="transparent", colour=NA),
          plot.background = element_rect(fill="transparent", color=NA), 
          panel.background = element_rect(fill="transparent", colour=NA))
  
  if(fc_line){
    p <- p +
      geom_vline(xintercept = fc_cutoff, 
                 linetype = "dashed",
                 color = line_colours[1], 
                 size = 0.5) +
      geom_vline(xintercept = -fc_cutoff, 
                 linetype = "dashed",
                 color = line_colours[1], 
                 size = 0.5)
  }
  
  # Add sig line if any significance
  if(any(toptable$padj <=  p_cutoff) & p_line ==  TRUE) {
    top <- max(toptable$pvalue[toptable$padj <=  p_cutoff], na.rm = TRUE)
    p <- p + geom_hline(yintercept = -log10(top), 
                        linetype = "dashed",
                        color = line_colours[2], 
                        size = 0.5)
  }
  
  # Add labels if any selected
  if(! is.null(label_rows)){
    if(! all(is.numeric(label_rows))) {
      if(! all(label_rows %in% rownames(toptable))){
        stop("label_rows must be in rownames(polar)")
      }}
    if(all(is.numeric(label_rows))) {
      if(! all(label_rows < nrow(toptable))){
        stop("label_rows not in 1:nrow(polar)")
      }}
    
    
    label_df <- toptable[label_rows, ]
    label_df$lp <- -log10(label_df$pvalue)
    
    p <- p + geom_text_repel(data = label_df, 
                             aes_string(x = "logFC", 
                                        y = "lp",
                                        label = "label"), 
                             size = label_size,
                             box.padding = unit(0.35, "lines"),
                             point.padding = unit(0.3, "lines"))
  }
  
  return(p)
}

#' Volcano Plots for a three-way comparison
#'
#' This function creates a volcano plot for all combinations of groups in a
#' factor.
#' @param polar A polar object with the pvalues between groups of interest. 
#' Created by \code{\link{polar_coords}}.
#' @param p_cutoff The cut-off for pvalue significance (default = 0.05). 
#' @param fc_cutoff The cut-off for fold change significance (default = 1). 
#' @param label_rows Row numbers or names for values to be annotated/labelled
#' (default = NULL).
#' @param label_size The font size of labels (default = 3)
#' @param text_size The font size of text (default = 10)
#' @param marker_size The size of markers (default = 0.7)
#' @param shared_legend_size The size for the legend (default = 1). 
#' @param sig_names A character vector of labels to be used for: 
#' non-significant; adjusted p < p_cutoff; |Fold Change| > fc_cutoff; and 
#' finally adjusted p < p_cutoff. 
#' If NULL c('Not Significant', \code{paste('Padj <', `p_cutoff`)},
#' \code{paste('|FC| >', `fc_cutoff`)},
#' \code{paste('Padj <', `p_cutoff`, 'and |FC| >', `fc_cutoff`)}) is used.
#' @param colours A character vector of colours to be used for non-significant;
#' adjusted p < p_cutoff;
#' |Fold Change| >, fc_cutoff; and adjusted p < p_cutoff. default = p_cutoff &
#' |Fold Change| > fc_cutoff markers respectively
#' (default = c('grey60', 'salmon', 'steelblue', 'limegreen')).
#' @param fc_line Logical whether to add vertical dashed line at fc_cutoff
#' (default = TRUE).
#' @param p_line Logical whether to add horizontal dashed line at p_cutoff
#' (default = TRUE).
#' @param line_colours A character vector stating the colour of lines to be
#' used for fc_line and p_line respectively (default = c('black', 'black')).
#' @param share_axes Logical whether plots should share axes when plotted 
#' together.
#' @importFrom ggplot2 layer_scales lims
#' @importFrom ggpubr ggarrange rremove get_legend as_ggplot
#' @importFrom stats p.adjust setNames
#' @return Returns a list of ggplot volcano plots. The first three elements 
#' contain comparisons between all contrasts. The last element in the list is a 
#' combined figure for all three plots.
#' @concept volcanoplot 
#' @keywords hplot
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @export
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(sampledata = syn_example_meta,
#'                           contrast = "Pathotype", 
#'                           groups = NULL, 
#'                           pvalues = syn_example_p, 
#'                           expression = syn_example_rld, 
#'                           p_col_suffix = "pvalue", 
#'                           padj_col_suffix = "padj", 
#'                           non_sig_name = "Not Significant", 
#'                           multi_group_prefix = "LRT",
#'                           significance_cutoff = 0.01, 
#'                           fc_col_suffix='log2FoldChange',
#'                           fc_cutoff = 0.3)
#' syn_volcano_plots <- volcano_trio(polar=syn_polar)
#' syn_volcano_plots$All

volcano_trio <- function(polar,
                         p_cutoff = 0.05,
                         fc_cutoff = 1,
                         label_rows = NULL,
                         label_size = 3,
                         text_size = 10,
                         marker_size = 0.7,
                         shared_legend_size = 1, 
                         sig_names = NULL,
                         colours = c("grey60", "salmon", 
                                     "steelblue", "limegreen"),
                         fc_line = TRUE,
                         p_line = TRUE,
                         line_colours = c("black", "black"),
                         share_axes = TRUE) {
  
  if(class(polar) != "polar") stop('polar must be an object of class polar')
  pvalues <- polar@pvalues
  
  if(length(colours) != 4) stop('The colour vector must be of length 4')
  if(any(unlist(lapply(colours, function(x) {
    class(try(col2rgb(x), silent = TRUE)) == "try-error"
  })))) {
    stop(paste(paste(colours[unlist(lapply(colours, function(x) {
      class(try(col2rgb(x), silent = TRUE)) == "try-error"
    }))], collapse=", "), 'is not a valid colour'))
  }
  
  if(length(line_colours) != 2) {
    stop('The line_colour vector must be of length 2 for fc and p cutoff
           lines respectively')
  }
  if( any(unlist(lapply(line_colours, function(x) {
    class(try(col2rgb(x), silent = TRUE)) == "try-error"
  }))))  {
    stop(paste(paste(line_colours[unlist(lapply(line_colours, function(x) {
      class(try(col2rgb(x), silent = TRUE)) == "try-error"
    }))], collapse=", "), 'is not a valid colour'))
  }
  
  if(! is.logical(fc_line )) {
    stop('fc_line must be a logical dictating whether to add vertical lines
           for fold change cutoff')
  }
  if(! is.logical(p_line )) {
    stop('p_line must be a logical dictating whether to add horizontal lines
           for p-value cutoff')
  }
  if(! is.logical(share_axes)) {
    stop('share_axes must be a logical dictating whether combined plots should 
         share axes')
  }
  
  v <- c("p_cutoff", "label_size", "fc_cutoff", "text_size", "marker_size", 
         "shared_legend_size")
  # Error if not numeric
  if( any(unlist(lapply(v, function(x) ! is.numeric(get(x))))) ) {
    stop(paste(
      paste(v[unlist(lapply(v, function(x) ! is.numeric(get(x))))], 
            collapse=", "),
      'must be numeric'))
  }
  
  # Error if out of ranges
  if(! (0 <= p_cutoff & p_cutoff <= 1)) stop('p_cutoff must be between 0 and 1')
  v <- c("label_size", "fc_cutoff", "text_size", "marker_size", 
         "shared_legend_size")
  if(any(unlist(lapply(v, function(x) get(x) <= 0)))) {
    stop(paste(v[unlist(lapply(v, function(x) get(x) <= 0))], 
               'must be greater than 0'))
  }
  
  # Extract the groups/levels
  groups <- gsub("_pvalue", "",
                 colnames(pvalues)[grepl("pvalue", colnames(pvalues))])
  if(! is.null(polar@multi_group_test)){ 
    groups <- groups[groups != polar@multi_group_test]
  }
  
  if(! all(paste0(groups, "_logFC") %in% colnames(pvalues)) ) {
    stop(paste('No', 
               paste(groups[ ! paste0(groups, "_logFC") %in% 
                               colnames(pvalues)], 
                     collapse=", "), 'fold-change columns'))
  }
  
  
  plot_outputs <- lapply(groups, function(x) {
    volcano_plot(pvalues_df = pvalues, 
                 comparison = x, 
                 p_cutoff = p_cutoff,
                 fc_cutoff = fc_cutoff,
                 label_rows = label_rows,
                 label_size = label_size,
                 text_size = text_size,
                 marker_size = marker_size,
                 shared_legend_size = shared_legend_size, 
                 sig_names = sig_names,
                 colours = colours,
                 fc_line = fc_line,
                 p_line = p_line,
                 line_colours = line_colours)
  })
  names(plot_outputs) <- groups
  
  # Output the plots using ggpubr and ggarrange
  # Output only one legend when plotting all
  out <- plot_outputs
  out[[1]] <- out[[1]] + rremove("legend")
  out[[2]] <- out[[2]] + rremove("legend")
  
  legend <- get_legend(out[[3]])
  
  out[[3]] <- out[[3]] + rremove("legend")
  out[[4]] <- as_ggplot(legend)
  
  if(share_axes){
    yRange <- unlist(lapply(out, function(x) layer_scales(x)$y$range$range))
    xRange <- unlist(lapply(out, function(x) layer_scales(x)$x$range$range))
    yMax <- max(yRange, na.rm = TRUE)
    xMax <- max(abs(xRange), na.rm = TRUE)
    
    out[1:3] <- lapply(out[1:3], function(p) {
      p + lims(x = c(-1*xMax, xMax), y = c(0, yMax))
    })
  }
  
  plot_outputs[["All"]] <- ggarrange(plotlist = out, 
                                     ncol = 4, 
                                     nrow = 1,
                                     widths = c(2, 2, 2, shared_legend_size))
  
  return(plot_outputs)
}

