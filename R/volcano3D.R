#' Three-Dimensional Volcano Plot
#'
#' Plots the pvalues from three-way comparisons in 3D space using plotly.
#' @param polar A polar object with created by \code{\link{polar_coords}}.
#' @param grid Optional grid list output by \code{\link{polar_grid}}. If 
#' NULL this will be calculated. 
#' @param fc_or_zscore whether to use fold change or z-score for the p-values 
#' (default = 'zscore').
#' @param label_rows character vector of rows (names or indices) to label.
#' @param label_size font size for labels (default = 14).
#' @param colours The named vector of colours for the groups. If NULL colours
#' will be assigned as c("green3", "cyan", "gold2", "blue", "purple" "red")
#' @param non_sig_colour The colour for non-significant markers according to 
#' fold change.
#' @param colour_scale whether to use a "discrete" or "continuous" colour scale 
#' (default = "discrete").
#' @param continuous_shift the number of degress (between 0 and 360) 
#' corresponding to the angle to offset the continuous colour scale by. The 
#' continuous colour scale is calculated by converting the angle to hue where 0 
#' degrees corresponds to red and 360 degrees to magenta (defualt = 120). 
#' @param axis_angle Angle in radians for the z axis (default = 0.5). 
#' @param z_aspectratio The aspect ratio for the z axis compared to x and y 
#' (default = 1). Decreasing this makes the plot appear more squat. 
#' @param xy_aspectratio The aspect ratio for the xy axis compared to z
#' (default = 1). Decreasing this makes the grid wider in the plot window. 
#' @param ... Optional parameters to pass to \code{\link[base]{polar_drid}} or 
#' \code{\link[base]{pretty}} to define the radial axis ticks.
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
#' @keywords pvalue, polar, plot
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
#' volcano3D(syn_polar, 
#'        label_rows = c("FMOD", "LAMP5", "TNNT3"), 
#'        label_size = 10, 
#'        xy_aspectratio = 1, 
#'        z_aspectratio = 0.9)


volcano3D <- function(polar,
                      fc_or_zscore = "zscore",
                      colours=NULL, 
                      non_sig_colour = "grey",
                      grid = NULL, 
                      label_rows = c(),
                      label_size = 14,
                      colour_scale = "discrete",
                      continuous_shift = 120, 
                      axis_angle = 0.5, 
                      z_aspectratio = 1, 
                      xy_aspectratio = 1, 
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
        volcano_toptable[, paste(polar@multi_group_test, "pvalue")])
    
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
    }
    
    polar_grid <- grid$polar_grid
    axes <- grid$axes
    axis_labels <- grid$axis_labs
    
    h <- as.numeric(grid$z)
    R <- as.numeric(grid$r)
    
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
    volcano_toptable$col[volcano_toptable$sig == polar@non_sig_name] <- 
        non_sig_colour
    
    cols <- setNames(as.character(unique(factor(volcano_toptable$col))), 
                     as.character(unique(droplevels(volcano_toptable$sig))))
    cols <- cols[match(levels(droplevels(polar_df$sig)), names(cols))]
    
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
            row <- volcano_toptable[volcano_toptable$Name == i, ]
            annot <- list(x = row$x, y = row$y, z = row$logP, text = row$Name, 
                          textangle = 0, ax = 75, ay = 0, arrowcolor = "black", 
                          arrowwidth = 1, arrowhead = 6, arrowsize = 1.5, 
                          xanchor = "left", yanchor = "middle")
        })
    } else {annot <- list()}
    
    axis_settings_xy[['range']] <- c(-1*(grid$r+1), grid$r+1)
    
    plot_ly(data = volcano_toptable, x = ~x, y = ~y, z = ~logP,
            marker = list(size = 2.6), 
            color = ~switch(colour_scale,
                            "discrete" = sig,
                            "continuous" = I(hue)),
            hoverinfo = 'text', 
            colors = switch(colour_scale,
                            "discrete" = cols,
                            "continuous" = NULL),
            text = ~paste0(rownames(volcano_toptable), 
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
        add_text(x = c(rep(1.05*R*sinpi(axis_angle), grid$n_z_breaks), 
                       1.2*R*sinpi(axis_angle)),
                 y = c(rep(1.05*R*cospi(axis_angle), grid$n_z_breaks), 
                       1.2*R*cospi(axis_angle)),
                 z = c(grid$z_breaks, h/2)*0.95,
                 text = c(grid$z_breaks, '-log<sub>10</sub>P'),
                 textposition = 'middle left', textfont = list(size = 10),  
                 color = I("black"), hoverinfo = 'none', showlegend = FALSE, 
                 inherit = FALSE) %>%
        
        add_trace(x = R*sinpi(axis_angle), y = R*cospi(axis_angle), 
                  z = c(0, h), color = I("black"),line = list(width = 2), 
                  showlegend = FALSE, type = "scatter3d", mode = "lines",
                  hoverinfo = "none", inherit = FALSE) %>%
        
        # label radial axis
        add_text(x = grid$text_coords$x, y = grid$text_coords$y, z = 0.05,
                 text = grid$text_coords$text, textposition = 'top center', 
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
               margin = list(l = 0, r = 0, b = 0, t = 0),
               scene = list(xaxis = axis_settings_xy,
                            yaxis = axis_settings_xy,
                            zaxis = axis_settings,
                            aspectratio = list(x = xy_aspectratio, 
                                               y = xy_aspectratio, 
                                               z = z_aspectratio)))
}
