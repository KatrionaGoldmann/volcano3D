#' Three-way radial comparison Polar Plot (using plotly)
#'
#' This function creates a plot using plotly which maps the expression
#' onto a polar coordinates.
#' @param polar The coordinates for plotting onto a polar plot
#' @param colours The named vector of colours for the groups. If NULL colours
#' will be assigned as c("green3", "cyan", "gold2", "blue", "purple" "red")
#' @param non_sig_colour The colour for non-significant markers according to 
#' fold change.
#' @param label_rows A character vector of genes or row numbers to label
#' @param grid Optional grid list output by \code{\link{polar_grid}}. If NULL 
#' this will be calculated. 
#' @param fc_or_zscore Whether to use the z-score or fold change as magnitude
#' (options are c("zscore", "fc"))
#' @param label_size Font size of labels (default = 14)
#' @param axis_title_size Font size for axis titles (default = 16)
#' @param axis_label_size Font size for axis labels (default = 10)
#' @param arrow_length length of label arrow (default = 50)
#' @param axis_ticks a numerical vector of tick breaks for the radial axis. If
#' NULL this will be calculated using \code{\link[base]{pretty}}.
#' @param axis_angle Angle for the axis labels (default = 5/6).
#' @param ... Optional parameters to pass to \code{\link[base]{pretty}} to 
#' define the radial axis ticks.
#' @return Returns a polar plotly plot featuring variables on a tri-axis
#' radial graph
#' @importFrom plotly plot_ly add_trace add_text add_markers layout
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @keywords pvalues, polar, plot, plotly, radial
#' @export
#' @examples
#' library(volcano3Ddata)
#' data(syn_data)
#' syn_p_obj <- create_dep(sampledata = syn_metadata, 
#'                     contrast = "Pathotype", 
#'                     pvalues = syn_pvalues,
#'                     p_col_suffix = "pvalue", 
#'                     fc_col_suffix = "log2FoldChange",
#'                     multi_group_prefix = "LRT", 
#'                     expression = syn_rld)
#' syn_polar <- polar_coords(dep = syn_p_obj)
#' radial_plotly(polar=syn_polar, 
#'             label_rows=c("SLAMF6", "PARP16", "ITM2C"))

radial_plotly <- function(polar,
                          colours=NULL, 
                          non_sig_colour = "grey60",
                          label_rows = NULL,
                          grid = NULL,
                          fc_or_zscore = "zscore",
                          label_size = 14,
                          axis_title_size = 16,
                          axis_label_size = 10,
                          arrow_length = 50,
                          axis_ticks = NULL,
                          axis_angle = 5/6,
                          ...){
    
    if(! class(polar) %in% c("polar")) stop("polar must be a polar object")
    polar_df <- polar@polar
    
    if(class(try(col2rgb(non_sig_colour),silent = TRUE)) == "try-error") {
        stop('non_sig_colour must be a valid colour')
        
    }
    if(any(unlist(lapply(colours, function(x) {
        class(try(col2rgb(x), silent = TRUE)) == "try-error"
    })))) stop('all values in colours must be valid colours')
    
    sig_levels <- levels(polar_df$sig)[levels(polar_df$sig) != 
                                          polar@non_sig_name]
    if(is.null(colours)){
        colours <- setNames(c("green3", "cyan", "gold2", "blue", 
                             "purple", "red"), 
                           sig_levels)
    }
    if(( ! is.null(names(colours) )) & 
       length(sig_levels[! sig_levels %in% names(colours)]) != 0) {
        stop(paste('No colour for', 
                   paste(sig_levels[! sig_levels %in% names(colours)], 
                         collapse=", ")))
    } 
    
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
    
    if(is.null(grid)) grid <- polar_grid(r_vector = polar_df$r, 
                                         axis_ticks = NULL,
                                         axis_angle = axis_angle, 
                                         ...)
    
    
    polar_grid <- grid$polar_grid
    axes <- grid$axes
    axis_labs <- grid$axis_labs
    r <- grid$r
    text_coords <- grid$text_coords
    
    
    # Set up the colours - pick the most highly expressed group
    polar_df$col <- as.character(colours[match(polar_df$sig, names(colours))])
    polar_df$col[polar_df$sig == polar@non_sig_name] <- non_sig_colour
  
    colours <- c("ns"=non_sig_colour, colours)
    names(colours)[names(colours) == "ns"] <- polar@non_sig_name
    colour_levels <- colours  
    
    # Align the levels
    polar_df$sig <- factor(polar_df$sig, levels=names(colour_levels))
    polar_df$col <- factor(polar_df$col, levels=colour_levels)
    
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
                 text = i,
                 textangle = 0,
                 ax = sign(row$x)*arrow_length*grid$r*cos(theta),
                 ay  = -1*sign(row$x)*arrow_length*grid$r*sin(theta),
                 font = list(color = gsub("[[:digit:]]+", "", row$col),
                             size = label_size),
                 arrowcolor = gsub("[[:digit:]]+", "", row$col),
                 arrowwidth = 1,
                 arrowhead = 6,
                 arrowsize = 1.5)
        })
    } else {annot <- list()}
    
    # Plotly plot
    p <- plot_ly(data = polar_df, x = ~x, mode = "none", type = "scatter",
                 colors = levels(polar_df$col), source = "BOTH", 
                 showlegend = FALSE) %>%
        #add the grid
        add_trace(x = polar_grid$x, y = polar_grid$y, color = I("#CBCBCB"),
                  line = list(width = 1), showlegend = FALSE, type = "scatter",
                  mode = "lines", hoverinfo = "none") %>%
        #add the "horizontal" L,M,F axes
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
               xaxis = list(title = "", zeroline = FALSE, showline = FALSE,
                            showticklabels = FALSE, showgrid = FALSE,
                            scaleratio = 1, scaleanchor = "y", 
                            autoscale = TRUE),
               yaxis = list(title = "", range = c(-0.5, 0.5), zeroline = FALSE,
                            showline = FALSE,	showticklabels = FALSE,
                            showgrid = FALSE),
               plot_bgcolor = "rgba(0,0,0,0)", autosize = TRUE,
               annotations = annot) %>%
        # label radial axis
        add_text(x = text_coords$x, y =  -text_coords$y, 
                 text = text_coords$text, textposition = 'top center', 
                 textfont = list(size = axis_label_size), color = I("black"), 
                 hoverinfo = 'none', showlegend = FALSE, inherit = FALSE) %>%
        # add the markers
        add_markers(data = polar_df, x = ~x, y = ~y,
                    marker = list(size = 6, sizemode = 'diameter'),
                    hoverinfo = 'text', key = rownames(polar_df), 
                    inherit = TRUE, color = ~sig, type = "scatter", 
                    mode = "markers", text = rownames(polar_df), 
                    showlegend = TRUE)
    
    return(p)
    
}
