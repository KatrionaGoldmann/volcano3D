#' Three-way radial comparison Polar Plot (using 'plotly')
#'
#' This function creates an interactive 'plotly' object which maps differential
#' expression onto a polar coordinates.
#' @param polar A polar object with the pvalues between groups of interest and
#' polar coordinates. Created by \code{\link{polar_coords}}.
#' @param colours A vector of colour names or hex triplets for each of the
#' six groups. Default = c("green3", "cyan", "blue",
#' "purple", "red", "gold2"). Colours are assigned in order: group1+,
#' group1+group2+, group2+, group2+group3+, group3+, group1+group3+.
#' @param non_sig_colour The colour for non-significant markers
#' (default = "grey60").
#' @param colour_scale whether to use a 'discrete' or 'continuous' colour scale
#' (default = 'discrete').
#' @param continuous_shift The number of radians (between 0 and 6)
#' to offset the continuous colour scale by. This is calculated by converting
#' the angle to a hue using \code{\link[grDevices]{hsv}} where 0 corresponds to
#' the colour scale starting with red and 360 with magenta (default = 2).
#' @param label_rows A vector of row names or numbers to label.
#' @param arrow_length The length of label arrows (default = 50).
#' @param grid An optional grid object. If NULL this will be calculated using
#' the default values of \code{\link{polar_grid}}.
#' @param fc_or_zscore Whether to use the z-score or fold change as magnitude.
#' Options are 'zscore' (default) or 'fc'.
#' @param label_size Font size of labels/annotations (default = 14)
#' @param axis_title_size Font size for axis titles (default = 16)
#' @param axis_label_size Font size for axis labels (default = 10)
#' @param axis_ticks A numerical vector of radial axis tick breaks. If
#' NULL this will be calculated using \code{\link[base]{pretty}}.
#' @param axis_angle Angle in radians for the radial axis (default = 5/6).
#' @param plot_height Plot height in px (default=700).
#' @param plot_width Plot width in px (default=700).
#' @param ... Optional parameters to pass to
#' \code{\link[volcano3D]{polar_grid}}.
#' @return Returns a 'plotly' plot featuring variables on a tri-axis
#' radial graph
#' @importFrom plotly plot_ly add_trace add_text add_markers layout
#' @importFrom stats p.adjust setNames
#' @importFrom grDevices hsv
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @keywords hplot iplot
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
#'                           fc_cutoff = 0.3)
#'
#' radial_plotly(polar = syn_polar, label_rows = c("SLAMF6"))


radial_plotly <- function(polar,
                          colours = c("green3", "cyan", "blue",
                                      "purple", "red", "gold2"),
                          non_sig_colour = "grey60",
                          colour_scale = "discrete",
                          continuous_shift = 1.33,
                          label_rows = NULL,
                          arrow_length = 50,
                          grid = NULL,
                          fc_or_zscore = "zscore",
                          label_size = 14,
                          axis_title_size = 16,
                          axis_label_size = 10,
                          axis_ticks = NULL,
                          axis_angle = 5/6,
                          plot_height = 700,
                          plot_width = 700,
                          ...){

    if(! class(polar) %in% c("polar")) stop("polar must be a polar object")
    polar_df <- polar@polar

    if(! is.null(colours) & length(colours) != 6){
        stop(paste("colours must be a character vector of plotting colours of",
                   "length six. One colour for each significance group"))
    }

    groups <- levels(polar@sampledata[, polar@contrast])
    sig_groups <- c(
        paste0(groups[1], "+"),
        paste0(groups[1], "+", groups[2], "+"),
        paste0(groups[2], "+"),
        paste0(groups[2], "+", groups[3], "+"),
        paste0(groups[3], "+"),
        paste0(groups[1], "+", groups[3], "+")
    )

    # If no colours selected dafault to rgb
    if(is.null(colours)){
        colours <- c("green3", "cyan", "blue", "purple", "red", "gold2")
    }
    if(is.null(non_sig_colour)){
        stop('Please enter a valid non_sig_colour')
    }

    # check if hex or can be converted to hex
    colours <- unlist(lapply(c(non_sig_colour, colours), function(x) {
        if(! grepl("#", x) &
           class(try(col2rgb(x), silent = TRUE))[1] == "try-error") {
            stop(paste(x, 'is not a valid colour'))
        } else if (! grepl("#", x) ) {
            y <- col2rgb(x)[, 1]
            x <- rgb(y[1], y[2], y[3], maxColorValue=255)
        }
        return(x)
    }))
    colours <- setNames(colours, c(polar@non_sig_name, sig_groups))

    if(! class(polar_df) %in% c("data.frame")) {
        stop("polar_df must be a data frame")
    }
    if(! fc_or_zscore %in% c("zscore", "fc")) {
        stop("fc_or_zscore must be either 'zscore' or 'fc'")
    }
    if(! is.numeric(label_size)) stop('label_size must be a numeric')
    if(! is.numeric(axis_title_size)) stop('axis_title_size must be a numeric')
    if(! is.numeric(axis_label_size)) stop('axis_label_size must be a numeric')
    if(! is.numeric(arrow_length)) stop('arrow_length must be a numeric')

    polar_df$x <- polar_df[, paste0("x_", fc_or_zscore)]
    polar_df$y <- polar_df[, paste0("y_", fc_or_zscore)]
    polar_df$r <- polar_df[, paste0("r_", fc_or_zscore)]
    polar_df <- polar_df[! is.nan(polar_df$x), ]

    if(is.null(grid)) {
        grid <- polar_grid(r_vector = polar_df$r,
                           r_axis_ticks = NULL,
                           axis_angle = axis_angle,
                           ...)
    } else{
        if(class(grid) != "grid") stop('grid must be a grid object')
    }
    if(! is.numeric(continuous_shift)) {
        stop('continuous_shift must be numeric')
    }
    if(! (0 <= continuous_shift & continuous_shift <= 2) ) {
        stop('continuous_shift must be between 0 and 2')
    }

    grid@polar_grid <- grid@polar_grid[grid@polar_grid$area != "cylinder", ]
    polar_grid <- grid@polar_grid
    axes <- grid@axes
    axis_labs <- grid@axis_labs
    r <- grid@r
    text_coords <- grid@text_coords

    # Set up the colours - pick the most highly expressed group
    polar_df$col <- as.character(colours[match(polar_df$sig, names(colours))])

    # Calculate the continuous colours
    offset <- (polar_df$angle[!is.na(polar_df$angle)] + continuous_shift/2)
    offset[offset > 1] <- offset[offset > 1] - 1
    polar_df$hue <- hsv(offset, 1, 1)
    polar_df$hue[polar_df$sig == polar@non_sig_name] <- non_sig_colour

    if(any(duplicated(colours))){
        warning(paste("Some colours are repeated. These will be compressed",
                      "into one significance group"))

        colours <- setNames(
            unique(colours),
            unlist(lapply(unique(colours), function(x) {
                groups_involved <- names(colours)[colours == x]
                # Filter onto those included only, unless empty
                if(any(groups_involved %in% polar_df$sig)){
                    groups_involved <- groups_involved[groups_involved %in%
                                                           polar_df$sig]
                }
                sub(",\\s+([^,]+)$", " or\n\\1",
                    paste(groups_involved, collapse=", /n"))
            })))
        polar_df$sig <- factor(names(colours)[match(polar_df$col, colours)])
    }


    colour_levels <- colours
    colour_levels <- colour_levels[names(colour_levels) %in%
                                       c(levels(polar_df$sig),
                                         polar@non_sig_name)]

    # Align the levels
    polar_df$sig <- factor(polar_df$sig, levels=names(colour_levels))

    # Annotate gene labels
    if (length(label_rows) != 0) {
        if(! all(is.numeric(label_rows))) {
            if(! all(label_rows %in% rownames(polar_df))) {
                stop("label_rows must be in rownames(polar_df)")
            }}
        if(all(is.numeric(label_rows))) {
            if(! all(label_rows < nrow(polar_df))) {
                stop("label_rows not in 1:nrow(polar_df)")
            }}
        annot <- lapply(label_rows, function(i) {
            row  <- polar_df[i, ]
            theta <- atan(row$y/row$x)
            list(x = row$x,
                 y = row$y,
                 text = as.character(row$label),
                 textangle = 0,
                 ax = sign(row$x)*arrow_length*grid@r*cos(theta),
                 ay  = -1*sign(row$x)*arrow_length*grid@r*sin(theta),
                 font = list(color = row$col, size = label_size),
                 arrowcolor = row$col,
                 arrowwidth = 1,
                 arrowhead = 6,
                 arrowsize = 1.5)
        })
    } else {annot <- list()}

    polar_df <- polar_df[c(which(polar_df$hue == non_sig_colour),
                           which(polar_df$hue != non_sig_colour)), ]

    # 'plotly' plot
    p <- plot_ly(data = polar_df, x = ~x, mode = "none", type = "scatter",
                 colors = switch(colour_scale,
                                 "discrete" = colour_levels,
                                 "continuous" = NULL),
                 source = "BOTH",
                 showlegend = FALSE) %>%
        # add the grid
        add_trace(x = polar_grid$x, y = polar_grid$y, color = I("#CBCBCB"),
                  line = list(width = 1), showlegend = FALSE, type = "scatter",
                  mode = "lines", hoverinfo = "none") %>%
        # add the "horizontal" axes
        add_trace(x = axes$x, y = axes$y, color = I("black"),
                  line = list(width = 2), showlegend = FALSE, type = "scatter",
                  mode = "lines", hoverinfo = "none", inherit = FALSE) %>%
        # add the label text
        add_text(x = axis_labs$x, y = axis_labs$y,
                 text = levels(polar@sampledata[, polar@contrast]),
                 color = I("black"), type = "scatter", mode = "text",
                 textfont = list(size = axis_title_size),
                 textposition = axis_labs$adjust, hoverinfo = 'none',
                 showlegend = FALSE, inherit = FALSE) %>%
        layout(showlegend = TRUE,
               xaxis = list(title = "", range = c(-grid@r, NULL),
                            zeroline = FALSE, showline = FALSE,
                            showticklabels = FALSE, showgrid = FALSE,
                            scaleratio = 1, scaleanchor = "y",
                            autoscale = TRUE),
               yaxis = list(title = "", range = c(-grid@r, grid@r),
                            zeroline = FALSE,
                            showline = FALSE,	showticklabels = FALSE,
                            showgrid = FALSE),
               plot_bgcolor = "rgba(0,0,0,0)",
               paper_bgcolor = 'rgba(0,0,0,0)',
               autosize = TRUE,
               annotations = annot)  %>%
        # label radial axis
        add_text(x = text_coords$x, y =  -text_coords$y,
                 text = text_coords$text, textposition = 'top center',
                 textfont = list(size = axis_label_size), color = I("black"),
                 hoverinfo = 'none', showlegend = FALSE, inherit = FALSE) %>%
        # add the markers
        add_markers(data = polar_df, x = ~x, y = ~y, key = ~label,
                    color = ~switch(colour_scale,
                                    "discrete" = sig,
                                    "continuous" = I(hue)),
                    colors = switch(colour_scale,
                                    "discrete" = levels(polar_df$col),
                                    "continuous" = NULL),
                    marker = list(size = 6, sizemode = 'diameter'),
                    hoverinfo = 'text', key = rownames(polar_df),
                    inherit = TRUE, type = "scatter",
                    mode = "markers", text = rownames(polar_df),
                    showlegend = (colour_scale == "discrete"))

    return(p)

}
