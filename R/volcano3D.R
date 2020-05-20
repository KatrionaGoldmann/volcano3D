#' Three-Dimensional Volcano Plot
#'
#' Plots the pvalues from three-way comparisons in 3D space using plotly.
#' @param polar A polar object with created by \code{\link{polar_coords}}.
#' @param colours A vector of colour names or hex triplets for each of the 
#' six groups. Default = c("green3", "cyan", "blue", 
#' "purple", "red", "gold2"). Colours are assigned in order: group1+, 
#' group1+group2+, group2+, group2+group3+, group3+, group1+group3+. 
#' @param non_sig_colour The colour for non-significant markers according to 
#' fold change (default='grey60').
#' @param colour_scale whether to use a 'discrete' or 'continuous' colour scale 
#' (default = 'discrete').
#' @param continuous_shift the number of degrees (between 0 and 360) 
#' corresponding to the angle to offset the continuous colour scale by. The 
#' continuous colour scale is calculated by converting the angle to hue where 
#' @param label_rows A vector of row names or numbers to label.
#' @param grid An optional grid object. If NULL this will be calculated using 
#' default values of  \code{\link{polar_grid}}. 
#' @param fc_or_zscore whether to use fold change or z-score for the p-values 
#' (default = 'zscore').
#' @param label_size font size for labels (default = 14).
#' @param axis_angle Angle in radians for the z axis (default = 0.5). 
#' @param z_aspectratio The aspect ratio for the z axis compared to x and y 
#' (default = 1). Decreasing this makes the plot appear more squat. 
#' @param xy_aspectratio The aspect ratio for the xy axis compared to z
#' (default = 1). Decreasing this makes the grid wider in the plot window. 
#' @param plot_height The plot height in px. Default=700. 
#' @param ... Optional parameters to pass to \code{\link{polar_grid}}.
#' @return Returns a cylindrical 3D plotly plot featuring variables on a 
#' tri-axis radial graph with the -log10(multi-group test p-value) on the 
#' z-axis
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and 
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @importFrom plotly plot_ly add_trace add_text layout %>%
#' @importFrom grDevices hsv rgb
#' @keywords hplot iplot
#' @concept volcanoplot 
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
#' 
#' volcano3D(syn_polar, 
#'     label_rows = c("FMOD", "LAMP5", "TNNT3"), 
#'     xy_aspectratio = 1, 
#'     label_size = 10, 
#'     z_aspectratio = 0.9)



volcano3D <- function(polar,
                      colours=c("green3", "cyan", "blue", 
                                "purple", "red", "gold2"), 
                      non_sig_colour = "grey60",
                      colour_scale = "discrete",
                      continuous_shift = 120, 
                      label_rows = c(),
                      grid = NULL, 
                      fc_or_zscore = "zscore",
                      label_size = 14,
                      axis_angle = 0.5, 
                      z_aspectratio = 1, 
                      xy_aspectratio = 1,
                      plot_height = 700,
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
    if(! colour_scale %in% c("discrete", "continuous")) {
        stop("colour_scale must be either 'discrete' or 'continuous'")
    }
    
    if(! is.numeric(label_size)) stop('label_size must be a numeric')
    if(! is.numeric(continuous_shift)) {
        stop('continuous_shift must be a numeric')
    }
    if(! (0 <= continuous_shift & continuous_shift <= 360) ) {
        stop('continuous_shift must be between 0 and 360')
    }
    
    # Calculate the colours by significance
    offset <- (polar_df$angle_degrees[!is.na(polar_df$angle_degrees)] + 
                   continuous_shift)/360
    offset[offset > 1] <- offset[offset > 1] - 1
    polar_df$hue <- hsv(offset, 1, 1)
    polar_df$hue[polar_df$sig == polar@non_sig_name] <- non_sig_colour
    
    volcano_toptable <- cbind(polar_df,
                              polar@pvalues[match(rownames(polar_df), 
                                                  rownames(polar@pvalues)), ])
    volcano_toptable$logP <- -1*log10(
        volcano_toptable[, paste0(polar@multi_group_test, "_pvalue")])
    
    volcano_toptable$x <- volcano_toptable[, paste0("x_", fc_or_zscore)]
    volcano_toptable$y <- volcano_toptable[, paste0("y_", fc_or_zscore)]
    volcano_toptable$r <- volcano_toptable[, paste0("r_", fc_or_zscore)]
    
    # Create the cylindrical grid for the volcano to sit in
    if(is.null(grid)){
        grid <- polar_grid(r_vector = volcano_toptable$r, 
                           z_vector = volcano_toptable$logP,
                           r_axis_ticks = NULL, 
                           z_axis_ticks = NULL, 
                           ...)
    } else { if(class(grid) != "grid") stop('grid must be a grid object')}
    
    polar_grid <- grid@polar_grid
    axes <- grid@axes
    axis_labels <- grid@axis_labs
    
    h <- as.numeric(grid@z)
    R <- as.numeric(grid@r)
    
    axis_settings <- list(title = "", zeroline = FALSE, showline = FALSE, 
                          showticklabels = FALSE, showgrid = FALSE, 
                          autotick = FALSE)
    axis_settings_xy <- list(title = "", zeroline = FALSE, showline = FALSE, 
                             showticklabels = FALSE, showgrid = FALSE, 
                             autotick = FALSE, showspikes = FALSE)
    
    volcano_toptable <- volcano_toptable[! is.na(volcano_toptable$x), ]
    volcano_toptable <- volcano_toptable[! is.na(volcano_toptable$y), ]
    volcano_toptable <- volcano_toptable[! is.na(volcano_toptable$logP), ]
    
    volcano_toptable$col <- as.character(colours[match(volcano_toptable$sig, 
                                                       names(colours))])
    if(any(duplicated(colours))){
        warning(paste("Some colours are repeated. These will be compressed", 
                      "into one significance group."))
        colours <- setNames(
            unique(colours), 
            unlist(lapply(unique(colours), function(x) {
                groups_involved <- names(colours)[colours == x]
                
                # Filter onto those included only, unless empty
                if(any(groups_involved %in% volcano_toptable$sig)){
                    groups_involved <- groups_involved[groups_involved %in% 
                                                           volcano_toptable$sig]
                }
                sub(",\\s+([^,]+)$", " or\n\\1", 
                    paste(groups_involved, collapse=", \n"))
            })))
        volcano_toptable$sig <- factor(names(colours)[match(
            volcano_toptable$col, colours)])
    }
    
    
    if(length(label_rows) != 0){
        
        if(! all(is.numeric(label_rows))) {
            if(! all(label_rows %in% rownames(volcano_toptable))) {
                stop("label_rows must be in rownames(polar_df)")
            }}
        if(all(is.numeric(label_rows))) {
            if(! all(label_rows < nrow(volcano_toptable))) {
                stop("label_rows not in 1:nrow(polar_df)")
            }}
        annot <- lapply(label_rows, function(i) {
            row <- volcano_toptable[i, ]
            annot <- list(x = row$x, y = row$y, z = row$logP, text = row$label, 
                          textangle = 0, ax = 75, ay = 0, arrowcolor = "black", 
                          arrowwidth = 1, arrowhead = 6, arrowsize = 1.5, 
                          xanchor = "left", yanchor = "middle")
        })
    } else {annot <- list()}
    
    axis_settings_xy[['range']] <- c(-1*(grid@r+1), grid@r+1)
    
    plot_ly(data = volcano_toptable, x = ~x, y = ~y, z = ~logP,
            marker = list(size = 2.6), 
            height = plot_height,
            key=~label,
            color = ~switch(colour_scale,
                            "discrete" = sig,
                            "continuous" = I(hue)),
            hoverinfo = 'text', 
            colors = switch(colour_scale,
                            "discrete" = colours,
                            "continuous" = NULL),
            text = ~paste0(label, 
                           "<br>Theta = ", as.integer(angle_degrees), "\u00B0",
                           "<br>r = ", formatC(r, digits = 3), 
                           "<br>P = ", format(logP, digits = 3, 
                                              scientific = 3)),
            type = "scatter3d", mode = "markers") %>%
        
        add_trace(x = polar_grid$x, y = polar_grid$y, z = polar_grid$z, 
                  color = I("#CBCBCB"), line = list(width = 2),
                  showlegend = FALSE, type = "scatter3d", mode = "lines", 
                  hoverinfo = "none",inherit = FALSE) %>%
        
        # Horizontal axes
        add_trace(x = axes$x, y = axes$y, z = 0, color = I("black"),
                  line = list(width = 2), showlegend = FALSE, 
                  type = "scatter3d", mode = "lines", hoverinfo = "none", 
                  inherit = FALSE) %>%
        add_text(x = axis_labels$x, y = axis_labels$y, z = 0, 
                 text = levels(polar@sampledata[, polar@contrast]),
                 color = I("black"), type = "scatter3d", mode = "text", 
                 textfont = list(size = 16),textposition = 'middle center', 
                 hoverinfo = 'none', showlegend = FALSE, inherit = FALSE) %>%
        # label z axis
        add_text(x = c(rep(1.05*R*sinpi(axis_angle), grid@n_z_breaks), 
                       1.2*R*sinpi(axis_angle)),
                 y = c(rep(1.05*R*cospi(axis_angle), grid@n_z_breaks), 
                       1.2*R*cospi(axis_angle)),
                 z = c(grid@z_breaks, h/2)*0.95,
                 text = c(grid@z_breaks, '-log<sub>10</sub>P'),
                 textposition = 'middle left', textfont = list(size = 10),  
                 color = I("black"), hoverinfo = 'none', showlegend = FALSE, 
                 inherit = FALSE) %>%
        
        add_trace(x = R*sinpi(axis_angle), y = R*cospi(axis_angle), 
                  z = c(0, h), color = I("black"),line = list(width = 2), 
                  showlegend = FALSE, type = "scatter3d", mode = "lines",
                  hoverinfo = "none", inherit = FALSE) %>%
        
        # label radial axis
        add_text(x = grid@text_coords$x, y = grid@text_coords$y, z = 0.05,
                 text = grid@text_coords$text, textposition = 'top center', 
                 textfont = list(size = 10), color = I("black"), 
                 hoverinfo = 'none', showlegend = FALSE, inherit = FALSE) %>%
        
        add_trace(type = 'scatter3d', mode = 'text', 
                  x = unlist(lapply(annot, "[[", "x")),
                  y = unlist(lapply(annot, "[[", "y")), 
                  z = unlist(lapply(annot, "[[", "z")),
                  text = unlist(lapply(annot, "[[", "text")), 
                  showlegend = FALSE, inherit = FALSE,
                  textfont = list(size = label_size))  %>%
        
        layout(showlegend = TRUE, dragmode = "turntable",
               margin = list(0, 0, 0, 0),
               paper_bgcolor = 'rgba(0, 0, 0, 0)',
               plot_bgcolor = 'rgba(0, 0, 0, 0)',
               scene = list(xaxis = axis_settings_xy,
                            yaxis = axis_settings_xy,
                            zaxis = axis_settings,
                            aspectratio = list(x = xy_aspectratio, 
                                               y = xy_aspectratio, 
                                               z = z_aspectratio)))
}
