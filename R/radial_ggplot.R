#' Ggplot for Three Way Polar Plot
#'
#' This function creates a radar plot using ggplot for a three-way comparison
#' @param polar The coordinates for plotting onto a polar plot.
#' @param label_rows A list of row names or indices in polar to be labelled
#' on the plot.
#' @param grid Optional grid list output by \code{\link{polar_grid}}. If NULL 
#' this will be calculated.
#' @param colours The named vector of colours for the groups. If NULL colours
#' will be assigned as c("green3", "cyan", "gold2", "blue", "purple", "red")
#' @param non_sig_colour The colour for non-significant markers according to 
#' fold change.
#' @param marker_alpha The alpha parameter for 
#' \code{\link[ggplot2]{geom_point}} (default = 0.7).
#' @param marker_size Size of the markers (default = 3).
#' @param colour_scale whether to use a "discrete" or "continuous" colour scale 
#' (default = "discrete").
#' @param continuous_shift the number of degrees (between 0 and 360) 
#' corresponding to the angle to offset the continuous colour scale by. The 
#' continuous colour scale is calculated by converting the angle to hue where 0 
#' degrees corresponds to red and 360 degrees to magenta (default = 120). 
#' @param label_size Font size of labels (default = 5).
#' @param axis_title_size Font size for axis titles (default = 5)
#' @param axis_label_size Font size for axis labels (default = 3)
#' @param fc_or_zscore Whether to use the z-score or fold change as magnitude
#' (options are c("zscore", "fc")).
#' @param axis_angle Angle for the axis labels (default = 1/6).
#' @param arrow_length length of label arrow (default = 1).
#' @param legend_size Size for the legend text (default = 20). 
#' @param ... Other optional parameters to pass to \code{\link{polar_grid}} or 
#' \code{\link[base]{pretty}} to define the radial axis ticks.
#' @return Returns a polar ggplot featuring variables on a tri-axis radial graph
#' @importFrom ggplot2 theme ggplot labs geom_path geom_path geom_text annotate 
#' geom_point scale_color_manual aes element_blank coord_fixed geom_segment
#' arrow unit
#' @importFrom graphics text 
#' @keywords pvalue, polar, plot, ggplot
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
#'                     p_col_suffix = "pvalue", 
#'                     fc_col_suffix = "log2FoldChange",
#'                     multi_group_prefix = "LRT", 
#'                     expression = syn_rld)
#' syn_polar <- polar_coords(dep = syn_p_obj)
#' 
#' radial_ggplot(polar = syn_polar, 
#'               colours=setNames(c("green3", "blue", "red", "gold2", 
#'               "purple", "cyan"),
#'               c("Fibroid", "Lymphoid", "Myeloid", "Fibroid+Myeloid+",  
#'               "Lymphoid+Myeloid+", "Fibroid+Lymphoid+")),
#'               marker_size = 1.5,
#'               label_size = 2.5,
#'               axis_label_size = 1.5, 
#'               axis_title_size = 2.5, 
#'               legend_size=5, 
#'               label_rows = c("SLAMF6", "PARP16", "ITM2C"))

radial_ggplot <- function(polar,
                          label_rows = NULL,
                          grid = NULL,
                          colours = NULL,
                          non_sig_colour = "grey60",
                          colour_scale = "discrete",
                          continuous_shift = 120, 
                          marker_alpha = 0.7,
                          marker_size = 3,
                          label_size = 5,
                          axis_title_size = 5,
                          axis_label_size = 3,
                          fc_or_zscore = "zscore",
                          axis_angle = 1/6,
                          arrow_length = 1,
                          legend_size = 20,
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
                              "purple", "red"), sig_levels)
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
    if(! class(polar) %in% c("polar")) stop("polar must be a polar object")
    if(! fc_or_zscore %in% c("zscore", "fc")) {
        stop("fc_or_zscore must be either 'zscore' or 'fc'")
    }
    
    if(! is.numeric(label_size)) stop('label_size must be a numeric')
    if(! is.numeric(axis_title_size)) stop('axis_title_size must be a numeric')
    if(! is.numeric(axis_label_size)) stop('axis_label_size must be a numeric')
    if(! is.numeric(arrow_length)) stop('arrow_length must be a numeric')
    if(! is.numeric(marker_size)) stop('marker_size must be a numeric')
    if(! is.numeric(marker_alpha)) stop('marker_alpha must be a numeric')
    if(! (marker_alpha >=  0 & marker_alpha <=  1)) {
        stop('marker_alpha must be between 0 and 1')
    }
    
    polar_df$x <- polar_df[, paste0("x_", fc_or_zscore)]
    polar_df$y <- polar_df[, paste0("y_", fc_or_zscore)]
    polar_df$r <- polar_df[, paste0("r_", fc_or_zscore)]
    
    # Set up the colours - pick the most highly expressed group
    polar_df$col <- as.character(colours[match(polar_df$sig, names(colours))])
    polar_df$col[polar_df$sig == polar@non_sig_name] <- non_sig_colour
    
    # Calculate the continuous colours
    offset <- (polar_df$angle_degrees[!is.na(polar_df$angle_degrees)] + 
                   continuous_shift)/360
    offset[offset > 1] <- offset[offset > 1] - 1
    polar_df$hue <- hsv(offset, 1, 1)
    polar_df$hue[polar_df$sig == polar@non_sig_name] <- non_sig_colour
    
    
    # make sure the non-sig markers are on the bottom - reshuffle the order
    polar_df$sig <- 
        factor(polar_df$sig,
               levels = c(polar@non_sig_name,
                          as.character(
                              unique(polar_df$sig[polar_df$sig !=  
                                                   polar@non_sig_name]))))
    colours <- c(colours, "ns"=non_sig_colour)
    names(colours)[names(colours) == "ns"] <- polar@non_sig_name
    cols <- colours[match(levels(droplevels(polar_df$sig)), names(colours))]
    
    if(is.null(grid)) grid <- polar_grid(r_vector = polar_df$r, 
                                         r_axis_ticks = NULL,
                                         axis_angle = axis_angle, 
                                         ...)
    
    grid$polar_grid <- grid$polar_grid[grid$polar_grid$area != "cylinder", ]
    
    
    if(! is.null(label_rows)){
        if(! all(is.numeric(label_rows))) {
            if(! all(label_rows %in% rownames(polar_df))){
                stop("label_rows must be in rownames(polar_df)")
            }}
        if(all(is.numeric(label_rows))) {
            if(! all(label_rows < nrow(polar_df))){
                stop("label_rows not in 1:nrow(polar_df)")
            }}
        annotation_df <- polar_df[label_rows, ]
        annotation_df$theta <- atan(annotation_df$y/annotation_df$x)
        annotation_df$xend  <- arrow_length*sign(annotation_df$x)*
            abs(grid$r*cos(annotation_df$theta))
        annotation_df$yend <- arrow_length*sign(annotation_df$y)*
            abs(grid$r*sin(annotation_df$theta))
    }
    
    # markers are plotted in order of rows so push ns to the bottom of plot
    polar_df <- polar_df[c(which(polar_df$sig ==  polar@non_sig_name),
                           which(polar_df$sig !=  polar@non_sig_name)), ]
    
    # alignment for text
    hadj <- -1*sign(grid$axis_labs$x)
    hadj[hadj ==  -1] <- 0
    
    p <- ggplot(polar_df, aes(x = polar_df$x, y = polar_df$y)) +
        labs(x = "", y = "", color = "") +
        
        # Concentric circles and radial spokes
        geom_path(data = grid$polar_grid, 
                  aes(x = grid$polar_grid$x, y = grid$polar_grid$y), 
                  alpha = 0.2) +
        
        # Three radial axes
        geom_path(data = grid$axes, aes(x = grid$axes$x, y = grid$axes$y)) +
        
        # radial axes ticks
        geom_text(data = grid$text_coords,
                  aes(x = grid$text_coords$x,
                      y = grid$text_coords$y,
                      label = grid$text_coords$text),
                  vjust = -1,
                  size = axis_label_size) +
        
        # Axes titles (three groups)
        annotate(geom = "text", x = grid$axis_labs$x,  y = grid$axis_labs$y,
                 hjust = hadj,
                 vjust = -1*sign(grid$axis_labs$y),
                 label = levels(polar@sampledata[, polar@contrast]),
                 color = "black", size = axis_title_size) +
        
        # Add markers
        geom_point(aes(color = switch(colour_scale, 
                                      "discrete"=polar_df$sig, 
                                      "continuous"=polar_df$hue)),
                   size = marker_size,
                   alpha = marker_alpha) +
        
        scale_color_manual(values = 
                               switch(colour_scale, 
                                      "discrete"=as.character(cols), 
                                      "continuous"=
                                          levels(factor(polar_df$hue)))) +
        
        # Set the background colour etc.
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.key = element_blank(),
              legend.justification = c(1, 1),
              legend.position = switch(colour_scale, 
                                       "discrete"=c(1, 1), 
                                       "continuous"="none"),
              legend.text = element_text(size = legend_size)) +
        
        # Fix the aspect ratio
        coord_fixed(ratio = 1,
                    xlim = c(-grid$r, grid$r*1.25),
                    ylim = c(-grid$r, grid$r))
    
    # Add any labelling desired
    if(! is.null(label_rows)){
        p <- p + geom_segment(data = annotation_df,
                              aes(x = annotation_df$x,
                                  y = annotation_df$y,
                                  xend = 0.9*annotation_df$xend,
                                  yend = 0.9*annotation_df$yend),
                              colour = annotation_df$hue, size = 0.5,
                              arrow = arrow(length = unit(0, "cm"))) +
            geom_text(data = annotation_df,
                      aes(x = 0.95*annotation_df$xend,
                          y = 0.95*annotation_df$yend,
                          label = rownames(annotation_df)),
                      color = annotation_df$hue, size = label_size) +
            geom_point(data = annotation_df,
                       aes(x = annotation_df$x, y = annotation_df$y),
                       shape = 1,
                       color = "black",
                       size = marker_size)
        
    }
    return(p)
}
