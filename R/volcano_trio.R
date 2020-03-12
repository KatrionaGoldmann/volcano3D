#' Volcanos Plots for a three-way comparison
#'
#' This function creates a volcano plot for all combinations of  groups in a
#' factor.
#' @param polar A polar object with the pvalues between groups of interest. 
#' Created by \code{\link{polar_coords}}.
#' @param p_cutoff The cut-off for pvalue significance (default = 0.05). 
#' @param fc_cutoff The cut-off for fold change significance (default = 1). 
#' @param label_rows Row IDs or rownames for values to be annotated/labelled
#' (default = NULL).
#' @param label_column The column ID or number to use as an annotation for
#' label_rows values (default = NULL). If left as NULL, and label_rows is 
#' selected, the rownames of pvalues will be used.
#' @param label_size The font size of labels (default = 3)
#' @param text_size The font size of text (default = 10)
#' @param marker_size The font size of markers (default = 0.7)
#' @param shared_legend_size Size for the legend (default = 1). 
#' @param sig_names A character vector of labels to be used for: 
#' non-significant; adjusted p < p_cutoff; |Fold Change| > fc_cutoff; and 
#' finally adjusted p < p_cutoff. 
#' default = p_cutoff & |Fold Change| > fc_cutoff markers 
#' respectively. 
#' (if NULL c('Not Significant', \code{paste('Padj <', `p_cutoff`)},
#' \code{paste('|FC| >', `fc_cutoff`)},
#' \code{paste('Padj <', `p_cutoff`, 'and |FC| >', `fc_cutoff`)}) is used).
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
#' @importFrom ggplot2 ggplot labs geom_path geom_path geom_text annotate 
#' geom_point scale_colour_manual theme aes theme_classic element_text 
#' geom_hline geom_vline unit layer_scales lims
#' @importFrom ggpubr ggarrange rremove get_legend as_ggplot
#' @importFrom ggrepel geom_text_repel
#' @return Returns a list of volcano plots. The first three elements contain 
#' comparisons between all combinations of the three levels in the comparison 
#' factor. The last element in the list is a combined figure for all three 
#' plots.
#' @examples
#' library(volcano3Ddata)
#' data(syn_data)
#' syn_p_obj <- polar_coords(sampledata = syn_metadata, 
#'                     contrast = "Pathotype", 
#'                     pvalues = syn_pvalues,
#'                     p_col_suffix="pvalue", 
#'                     fc_col_suffix = "log2FoldChange",
#'                     multi_group_prefix = "LRT", 
#'                     expression=syn_rld)
#' syn_mod_plots <- volcano_trio(polar=syn_p_obj)
#' 
#' syn_mod_plots$All
#' syn_mod_plots$`Lymphoid-Myeloid`
#' @keywords volcano, pvalues, plot
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @export


volcano_trio <- function(polar,
                         p_cutoff = 0.05,
                         fc_cutoff = 1,
                         label_rows = NULL,
                         label_column = NULL,
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

  if(! is.null(label_column)) if(! any(grepl(label_column, colnames(pvalues)))){
    stop('label_column must be a column name in pvlaues refering to the rows
           to label')
  }
  if(length(colours) != 4) stop('The colour vector must be of length 4')
  if(length(line_colours) != 2) {
    stop('The line_colour vector must be of length 2 for fc and p cutoff
           lines respectively')
  }
  if(! is.logical(fc_line )) {
    stop('fc_line must be a logical dictating whether to add vertical lines
           for fold change cutoff')
  }
  if(! is.logical(p_line )) {
    stop('p_line must be a logical dictating whether to add horizontal lines
           for p-value cutoff')
  }
  if(! is.numeric(p_cutoff)) stop('p_cutoff must be a numeric')
  if(! is.numeric(fc_cutoff)) stop('fc_cutoff must be a numeric')
  if(! is.numeric(label_size)) stop('label_size must be a numeric')
  
  
  # Extract the groups/levels
  groups <- gsub(" pvalue", "",
                 colnames(pvalues)[grepl("pvalue", colnames(pvalues))])
  if(! is.null(polar@multi_group_test)){ 
    groups <- groups[groups != polar@multi_group_test]
  }
  
  if(! all(paste(groups, "logFC") %in% colnames(pvalues)) ) {
    stop(paste('No', 
               paste(groups[ ! paste(groups, "logFC") %in% colnames(pvalues)], 
                     collapse=", "), 'fold-change columns'))
  }
  
  # Create a list of volcano plots showing DEG between groups
  volcano_plot <- function(comparison){
    toptable <- pvalues[, grepl(comparison, colnames(pvalues))]
    
    if(is.null(sig_names)){
      sig_names <- c("Not Significant",
                     paste("Padj <", p_cutoff),
                     paste("|FC| >", fc_cutoff),
                     paste("Padj <", p_cutoff, "and |FC| >", fc_cutoff))
    }
    
    colnames(toptable)[unlist(lapply(c("pvalue", "padj", "logFC"), function(x) {
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
    
    # Create the volcano plot ggplot
    p <- ggplot(toptable, 
                aes(x = toptable$logFC, 
                    y = -log10(toptable$pvalue)),
                color = toptable$cols) +
      geom_point(aes(color = toptable$cols), size=marker_size) +
      scale_colour_manual(values = names(mapping)[match(levels(toptable$cols),
                                                        mapping)]) +
      labs(y = expression(-log["10"]*p), 
           x = expression(log["2"]*FC),
           title = gsub("-", " vs ", comparison), 
           color = "Significance") +
      theme_classic() +
      theme(text = element_text(size = text_size))
    
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
      if(is.null(label_column)) {
        labs <- rownames(toptable[label_rows, ])
        label_df <- toptable[label_rows, ]
      } else {
        labs <- toptable[label_rows, label_column]
        label_df <- toptable[toptable[, label_column] %in% label_rows]
      }
      
      p <- p + geom_text_repel(data = label_df, 
                                        aes(x = label_df$logFC, 
                                            y = -log10(label_df$pvalue),
                                            label = labs), 
                                        size = label_size,
                                        box.padding = unit(0.35, "lines"),
                                        point.padding = unit(0.3, "lines"))
    }
    
    return(p)
  }
  plot_outputs <- lapply(groups, volcano_plot)
  
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

