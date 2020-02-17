#' Ggplot for Three Way Polar Plot
#'
#' This function creates a radar plot using ggplot for a three-way comparison
#' @param polar The coordinates for plotting onto a polar plot.
#' @param label_rows A list of row names or indices in polar to be labelled
#' on the plot.
#' @param grid Optional grid list output by \code{\link{polar_grid}}. If NULL 
#' this will be calculated.
#' @param fc_cutoff The cut-off for fold change, below which markers will be
#' coloured according to `non_sig_colour`` (default = 0.3).
#' @param non_sig_colour The colour for non-significant markers according to 
#' fold change.
#' @param marker_alpha The alpha parameter for 
#' \code{\link[ggplot2]{geom_point}} (default = 0.7).
#' @param marker_size Size of the markers (default = 3).
#' @param label_size Font size of labels (default = 5).
#' @param axis_title_size Font size for axis titles (default = 5)
#' @param axis_label_size Font size for axis labels (default = 3)
#' @param fc_or_zscore Whether to use the z-score or fold change as magnitude
#' (options are c("zscore", "fc")).
#' @param axis_angle Angle for the axis labels (default = 1/6).
#' @param arrow_length length of label arrow (default = 1).
#' @param legend_size Size for the legend text (default = 20). 
#' @param ... Other optional parameters to pass to \code{\link{polar_grid}}
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
#' data(syn_data)
#' syn_p_obj <- create_dep(sampledata = syn_metadata, 
#'                     contrast = "Pathotype", 
#'                     pvalues = syn_pvalues,
#'                     p_col_suffix = "pvalue", 
#'                     fc_col_suffix = "log2FoldChange",
#'                     multi_group_prefix = "LRT", 
#'                     expression = syn_rld)
#' syn_polar <- polar_coords(dep = syn_p_obj)
#' radial_ggplot(polar = syn_polar, 
#'               fc_cutoff = 0.1, 
#'               label_rows = c("SLAMF6", "PARP16", "ITM2C")) 

radial_ggplot <- function(polar,
                          label_rows = NULL,
                          grid = NULL,
                          fc_cutoff = 0.3,
                          non_sig_colour = "grey60",
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
    
    
    polar_df <- polar@polar
    
    if(! class(polar_df) %in% c("data.frame")) {
        stop("polar_df must be a data frame")
    }
    if(! class(polar) %in% c("polar")) stop("polar must be a polar object")
    if(! fc_or_zscore %in% c("zscore", "fc")) {
        stop("fc_or_zscore must be either 'zscore' or 'fc'")
    }
    if(! is.numeric(fc_cutoff)) stop('fc_cutoff must be a numeric')
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
    
    polar_df$col[polar_df$r < fc_cutoff] <- non_sig_colour
    polar_df$sig[polar_df$r < fc_cutoff] <- polar@non_sig_name
    
    # make sure the non-sig markers are on the bottom - reshuffle the order
    polar_df$sig <- factor(polar_df$sig,
                          levels = c(polar@non_sig_name,
                                     as.character(
                                         unique(
                                             polar_df$sig[polar_df$sig !=  
                                                    polar@non_sig_name]))))
    cols <- setNames(as.character(unique(droplevels(polar_df$col))),
                     as.character(unique(droplevels(polar_df$sig))))
    cols <- cols[match(levels(droplevels(polar_df$sig)), names(cols))]
    
    if(is.null(grid)) grid <- polar_grid(r_vector = polar_df$r, 
                                         axis_ticks = NULL,
                                         axis_angle = axis_angle, 
                                         ...)
    
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
        geom_path(data = grid$axes,
                  aes(x = grid$axes$x, y = grid$axes$y)) +
          
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
        geom_point(aes(color = polar_df$sig),
                   size = marker_size,
                   alpha = marker_alpha) +
        
        scale_color_manual(values = as.character(cols)) +

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
              legend.position = c(1, 1),
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
                             colour = annotation_df$col, size = 0.5,
                             arrow = arrow(length = unit(0, "cm"))) +
            geom_text(data = annotation_df,
                      aes(x = 0.95*annotation_df$xend,
                          y = 0.95*annotation_df$yend,
                          label = rownames(annotation_df)),
                      color = annotation_df$col, size = label_size) +
            geom_point(data = annotation_df,
                       aes(x = annotation_df$x, y = annotation_df$y),
                       shape = 1,
                       color = "black",
                       size = marker_size)

    }
    return(p)
}
