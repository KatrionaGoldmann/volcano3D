#' Boxplot to compare groups
#'
#' Plots the expression of a specific row in expression to compare the three
#' groups in a boxplot using either ggplot or plotly.
#' @param polar A 'volc3d' object including expression data from groups of
#' interest. Created by \code{\link{polar_coords}}.
#' @param value The column name or number in \code{polar@data} to be
#' analysed
#' @param box_colours The fill colours for each box assigned in order of
#' levels_order. Default = c('green3', 'blue', 'red') ).
#' @param test The statistical test used to compare expression.
#' Allowed values include: \itemize{
#'  \item \code{polar_pvalue} (default) and 'polar_padj' for the pvalues
#'  and adjusted pvalues in the polar object.
#'  \item \code{polar_multi_pvalue} and \code{polar_multi_padj} for the pvalues
#'  and adjusted pvalues across all groups using the
#'  \code{polar@multi_group_test } columns.
#'  \item \code{\link[stats]{t.test}} (parametric) and
#'  \code{\link[stats]{wilcox.test}} (non-parametric). Perform comparison
#'  between groups of samples.
#'  \item \code{\link[stats]{anova}} (parametric) and
#'  \code{\link[stats]{kruskal.test}} (non-parametric). Perform one-way ANOVA
#'  test comparing multiple groups. }
#' @param levels_order A character vector stating the contrast groups to be
#' plotted, in order. If `NULL` this defaults to the levels in
#' `polar@outcome`.
#' @param my_comparisons A list of contrasts to pass to
#' \code{\link[ggpubr]{stat_compare_means}}. If `NULL` (default) all contrast
#' pvalues are calculated and plotted.
#' @param text_size The font size of text (default = 10)
#' @param stat_colour Colour to print statistics (default="black").
#' @param stat_size The font size of statistical parameter (default = 3).
#' @param step_increase The distance between statistics on the y-axis
#' (default = 0.05).
#' @param plot_method Whether to use 'plotly' or 'ggplot'. Default is 'ggplot'
#' @param ... Other parameters for \code{\link[ggpubr]{stat_compare_means}}
#' @return Returns a boxplot featuring the differential expression
#' between groups in comparison with annotated pvalues.
#' @importFrom ggpubr compare_means ggboxplot stat_pvalue_manual
#' stat_compare_means
#' @importFrom plotly layout plot_ly add_trace add_markers
#' @importFrom utils combn
#' @importFrom grDevices hsv col2rgb
#' @importFrom ggplot2 theme ggplot labs geom_path geom_path geom_text annotate
#' geom_point scale_color_manual aes geom_jitter element_rect aes_string
#' element_text
#' @importFrom methods is
#' @keywords hplot
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://pubmed.ncbi.nlm.nih.gov/31461658/}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @export
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
#'                           data = t(syn_example_rld))
#'
#' boxplot_trio(syn_polar, value = "COBL", plot_method="plotly")
#' boxplot_trio(syn_polar, value = "COBL")

boxplot_trio <- function(polar,
                         value,
                         box_colours = c('green3', 'blue', 'red'),
                         test = "polar_pvalue",
                         levels_order = NULL,
                         my_comparisons = NULL,
                         text_size = 10,
                         stat_colour = "black",
                         stat_size = 3,
                         step_increase = 0.05,
                         plot_method="ggplot",
                         ...){
  if (is(polar, "polar")) {
    args <- as.list(match.call())[-1]
    return(do.call(boxplot_trio_v1, args))  # for back compatibility
  }
  if(! is(polar, "volc3d")) stop("Not a 'volc3d' class object")
  outcome <- polar@outcome
  expression <- t(polar@data)
  pvalues <- polar@pvals
  padj <- polar@padj

  if(! test %in% c("polar_pvalue", "polar_padj", "polar_multi_pvalue",
                   "polar_multi_padj", "t.test", "wilcox.test", "anova",
                   "kruskal.test")) {
    stop(paste("expression must be a data frame or c('polar_pvalues',",
               "'polar_padj', 'polar_multi_pvalue', 'polar_multi_padj',",
               "'t.test', 'wilcox.test', 'anova', 'kruskal.test')"))
  }
  if(is.null(levels_order)) {
    levels_order <- levels(outcome)
  }
  if(! all(levels_order %in% levels(outcome))){
    stop(paste('levels_order must be a character vector defining the order',
               "of levels in 'outcome'"))
  }
  if(length(box_colours) != length(levels_order)){
    stop(paste0('The length of box_colours must match teh length of',
                'levels_order'))
  }
  # Convert to hex colours for plotly
  box_colours <- unlist(lapply(box_colours, function(x) {
    if(! grepl("#", x) &
       class(try(col2rgb(x), silent = TRUE))[1] == "try-error") {
      stop(paste(x, 'is not a valid colour'))
    } else if (! grepl("#", x) ) {
      y <- col2rgb(x)[, 1]
      x <- rgb(y[1], y[2], y[3], maxColorValue=255)
    }
    return(x)
  }))

  if(length(value) > 1) stop("value must be of length 1")
  if(! value %in% rownames(expression)) {
    stop("value/gene is not in rownames(expression)")
  }

  if(! plot_method %in% c('plotly', 'ggplot')){
    stop("plot_method must be either plotly or ggplot")
  }

  colour_map <- setNames(box_colours, levels_order)

  if(is.character(value)) {
    index <- which(rownames(expression) ==  value)
  }

  if(is.null(my_comparisons)) {
    comps <- levels_order
    my_comparisons <- lapply(seq_len(ncol(combn(comps, 2))), function(i) {
      as.character(combn(comps, 2)[, i])
    })
  }

  df <- data.frame(
                   "group" = outcome,
                   "row" = expression[value, ])
  df <- df[! is.na(df$row), ]
  df <- df[df$group %in% levels_order, ]

  # relevel based on defined order
  df$group <- factor(df$group, levels_order)
  df$col <- factor(df$group, labels=colour_map[match(levels(df$group),
                                                     names(colour_map))])

  map_pos <- setNames(seq_along(levels(df$group)), levels(df$group))

  # Calculate the pvalues depending on test
  if(test %in% c("t.test", "wilcox.test", "anova", "kruskal.test")){
    pvals <- compare_means(formula = row ~ group, data = df,
                           comparisons = my_comparisons,
                           method = test,
                           step.increase = step_increase,
                           size=stat_size)
    if (test %in% c("t.test", "wilcox.test")) {
      pvals$x.position <- map_pos[pvals$group1] +
        (map_pos[pvals$group2] - map_pos[pvals$group1])/2
      pvals$y.position <- max(df$row, na.rm=TRUE)*
        (1.01 + step_increase*c(seq_len(nrow(pvals))-1))
      pvals$new_p_label <- pvals$p.format
    } else pvals <- pvals$p
    
    # groups comparisons
  } else if (! grepl("multi", test)){
    
    pvals <- switch(test,
                    "polar_pvalue" = pvalues[value, 2:4],
                    "polar_padj" = padj[value, 2:4])
    pvals <- data.frame(p = pvals)
    pvals$group1 <- levels(outcome)[c(1,1,2)]
    pvals$group2 <- levels(outcome)[c(2,3,3)]
    pvals$p.format <- format(pvals$p, digits=2)
    pvals$method <- test
    pvals$y.position <- max(df$row, na.rm=TRUE)
    pvals$comp <- paste0(pvals$group1, "_", pvals$group2)
    pvals$x.position <- map_pos[pvals$group1] +
      (map_pos[pvals$group2] - map_pos[pvals$group1])/2
    pvals$y.position <- pvals$y.position[1]*
      (1.01 + step_increase*(seq_len(nrow(pvals))-1))
    comp_use <- unlist(lapply(my_comparisons, function(x) { 
      c(paste0(x[1], "_", x[2]), paste0(x[2], "_", x[1]))
      }))
    pvals <- pvals[pvals$comp %in% comp_use, ]

    # multi group comparisons
  } else{
    # polar_multi_test
    pvals <- switch(test,
                    "polar_multi_pvalue" = pvalues[value, 1],
                    "polar_multi_padj" = padj[value, 1])
  }

  if(plot_method == 'ggplot'){
    p <- ggboxplot(data = df,
                   x = "group",
                   y = "row",
                   xlab = "",
                   ylab = rownames(polar@pvals)[index],
                   fill = "group",
                   color = "group",
                   palette = box_colours,
                   outlier.shape = NA,
                   alpha = 0.3) +
      geom_jitter(data=df, height = 0, width = 0.30,
                  aes_string(color="group")) +
      theme(legend.position = "none",
            text = element_text(size = text_size),
            plot.background = element_rect(fill="transparent", color=NA),
            panel.background = element_rect(fill="transparent", colour=NA),
            legend.background = element_rect(fill="transparent", colour=NA))


    if(! grepl("multi", test) & !(test %in% c("anova", "kruskal.test"))){

      p <- p + stat_pvalue_manual(
        data = pvals, label = "p.format",
        xmin = "group1", xmax = "group2",
        step.increase = step_increase,
        y.position = "y.position", color = stat_colour,
        size=stat_size, ...)
    } else{
      # multi group comparison
      p <- p + annotate("text", x = 0.5 + length(unique(df$group))/2,
                        y = Inf, vjust = 2, hjust = 0.5, color = stat_colour,
                        label = paste("p =", format(pvals, digits = 2)))
    }
  } else{
    p <- df %>%
      plot_ly() %>%
      add_trace(x = ~as.numeric(group),  y = ~row,
                type = "box",
                colors = levels(df$col), color = ~col,
                opacity=0.5, marker = list(opacity = 0),
                hoverinfo="none", showlegend = FALSE) %>%
      add_markers(x = ~jitter(as.numeric(group)), y = ~row,
                  marker = list(size = 6, color=~col),
                  hoverinfo = "text",
                  text = paste0(rownames(df),
                    "<br>Group: ", df$group,
                                 "<br>Expression: ", format(df$row, digits = 3)),
                  showlegend = FALSE) %>%
      layout(legend = list(orientation = "h",
                           x =0.5, xanchor = "center",
                           y = 1, yanchor = "bottom"
      ),
      xaxis = list(title = "", tickvals = 1:3,
                   ticktext = levels(df$group)),
      yaxis = list(title = value))


    lines <- list()
    if(! grepl("multi", test)){
      for (i in seq_len(nrow(pvals))) {
        line <- list(line=list(color = stat_colour))
        line[["x0"]] <- map_pos[pvals$group1][i]
        line[["x1"]] <- map_pos[pvals$group2][i]
        line[c("y0", "y1")] <- pvals$y.position[i]
        lines <- c(lines, list(line))
      }

      a <- list(
        x = as.numeric(pvals$x.position),
        y = pvals$y.position,
        text = format(pvals$p.format, digits=3),
        xref = "x",
        yref = "y",
        yanchor = "bottom",
        font = list(color = stat_colour),
        showarrow = FALSE
      )

      p <- p %>% layout(annotations = a, shapes=lines)

    } else{
      a <- list(
        x = 2,
        y = step_increase + max(df$row, na.rm=TRUE),
        text = format(pvals, digits=3),
        xref = "x",
        yref = "y",
        font = list(color = stat_colour),
        showarrow = FALSE
      )

      p <- p %>% layout(annotations = a)
    }


  }

  return(p)
}

