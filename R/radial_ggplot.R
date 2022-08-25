#' 'Ggplot' for Three Way Polar Plot
#'
#' This function creates a 3-way polar plot using 'ggplot' for a three-class
#' comparison.
#' 
#' @param polar A 'volc3d' object with the p-values between groups of interest
#'   and polar coordinates created by \code{\link{polar_coords}},
#'   \code{\link{deseq_polar}} or \code{\link{voom_polar}}.
#' @param type Numeric value whether to use scaled (z-score) or unscaled (fold
#'   change) as magnitude. Options are 1 = z-score (default) or 2 =
#'   unscaled/fold change.
#' @param colours A vector of colours for the non-significant points and each of
#'   the six groups.
#' @param label_rows A vector of row names or indices to label
#' @param arrow_length The length of label arrows
#' @param label_size Font size of labels/annotations (default = 5).
#' @param colour_code_labels Logical whether label annotations should be colour
#' coded. If FALSE `label_colour` is used.
#' @param label_colour Colour of annotation labels if not colour coded
#' @param grid_colour The colour of the grid (default="grey80")
#' @param grid_width The width of the axis lines (default=0.6)
#' @param axis_colour The colour of the grid axes and labels (default="black")
#' @param axis_width The width of the axis lines (default=1)
#' @param axis_title_size Font size for axis titles (default = 5)
#' @param axis_label_size Font size for axis labels (default = 3)
#' @param marker_alpha The alpha parameter for markers (default = 0.7)
#' @param marker_size Size of the markers (default = 3)
#' @param marker_outline_colour Colour for marker outline (default = white)
#' @param marker_outline_width Width for marker outline (default = 0.5)
#' @param legend_size Size for the legend text (default = 20).
#' @param ... Optional parameters passed to \code{\link[volcano3D]{polar_grid}}
#'   e.g. `r_axis_ticks` or `axis_angle`
#' @return Returns a polar 'ggplot' object featuring variables on a tri-axis
#' radial graph
#' @importFrom ggplot2 theme ggplot labs geom_path geom_path geom_text annotate
#' geom_point scale_color_manual aes element_blank coord_fixed geom_segment
#' arrow unit element_rect aes_string scale_fill_manual element_text
#' @importFrom methods is
#' @keywords hplot
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://pubmed.ncbi.nlm.nih.gov/31461658/}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @seealso \code{\link{polar_coords}}
#' @export
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
#'                           data = t(syn_example_rld))
#'
#' radial_ggplot(polar = syn_polar, label_rows = c("COBL"))

radial_ggplot <- function(polar,
                          type = 1,
                          colours = NULL,
                          label_rows = NULL,
                          arrow_length = 1,
                          label_size = 5,
                          colour_code_labels = FALSE,
                          label_colour = "black",
                          grid_colour = "grey80", 
                          grid_width = 0.7,
                          axis_colour = "black",
                          axis_width = 1,
                          axis_title_size = 5,
                          axis_label_size = 3,
                          marker_alpha = 0.7,
                          marker_size = 3,
                          marker_outline_colour = "white",
                          marker_outline_width = 0.5,
                          legend_size = 20,
                          ...){
    if (is(polar, "polar")) {
      args <- as.list(match.call())[-1]
      return(do.call(radial_ggplot_v1, args))  # for back compatibility
    }
    if(! is(polar, "volc3d")) stop("polar must be a 'volc3d' object")
    polar_df <- polar@df[[type]]

    grid <- polar_grid(r_vector = polar_df$r, ...)

    grid@polar_grid <- grid@polar_grid[grid@polar_grid$area != "cylinder", ]

    # markers are plotted in order of rows so push ns to the bottom of plot
    polar_df <- polar_df[order(polar_df$lab), ]
    
    old_levels <- levels(polar_df$lab)
    polar_df$lab <- droplevels(polar_df$lab)
    if (is.null(colours)) {
      colours <- polar@scheme
      colours <- colours[old_levels %in% levels(polar_df$lab)]
    }

    # alignment for text
    hadj <- -1*sign(grid@axis_labs$x)
    hadj[hadj ==  -1] <- 0
    
    polar_df$cg <- polar_df$col

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
        annotation_df$label <- rownames(annotation_df)
        annotation_df$theta <- atan(annotation_df$y/annotation_df$x)
        annotation_df$xend  <- arrow_length*sign(annotation_df$x)*
            abs(grid@r*cos(annotation_df$theta))
        annotation_df$yend <- arrow_length*sign(annotation_df$y)*
            abs(grid@r*sin(annotation_df$theta))
    }
    
    p <- ggplot(polar_df, aes_string(x = "x", y = "y")) +
        labs(x = "", y = "", color = "") 

    # Concentric circles and radial spokes
    rem <- which(is.na(grid@polar_grid$x))
    invisible(lapply(1:(length(rem)-1), function(g){
        p <<- p + geom_path(data = grid@polar_grid[(rem[g]+1):(rem[g+1]-1), ],
                  aes_string(x = "x", y = "y"),
                  size=grid_width,
                  colour=grid_colour)
        }))
    
    # Three radial axes
    rem <- c(0, which(is.na(grid@axes$x) ))
    invisible(lapply(1:(length(rem)-1), function(g){
        p <<- p + geom_path(data = grid@axes[(rem[g]+1):(rem[g+1]-1), ],
                            aes_string(x = "x", y = "y"), 
                            size=axis_width, color=axis_colour)
    }))
    

    p <- p +

        # radial axes ticks
        geom_text(data = grid@text_coords,
                  aes_string(x = "x",
                             y = "y",
                             label = "text"),
                  color = axis_colour,
                  vjust = -1,
                  size = axis_label_size) +

        # Axes titles (three groups)
        annotate(geom = "text",
                 x = grid@axis_labs$x,
                 y = grid@axis_labs$y,
                 hjust = hadj,
                 vjust = -1*sign(grid@axis_labs$y),
                 label = levels(polar@outcome),
                 color = axis_colour, size = axis_title_size) +

        # Add markers
        geom_point(aes_string(fill = "lab"),
                   size = marker_size,
                   alpha = marker_alpha, 
                   color = marker_outline_colour,
                   stroke = marker_outline_width,
                   shape = 21) +

        scale_fill_manual(name="",
            values = as.character(colours)) +
          
        # Set the background colour etc.
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.key = element_blank(),
              legend.justification = c(1, 1),
              legend.position = "right",
              legend.text = element_text(size = legend_size),
              legend.background = element_rect(fill="transparent", colour=NA),
              plot.background = element_rect(fill="transparent", color=NA),
              panel.background = element_rect(fill="transparent", colour=NA)) +

        # Fix the aspect ratio
        coord_fixed(ratio = 1,
                    xlim = c(-grid@r, grid@r*1.25),
                    ylim = c(-grid@r, grid@r))
    
    # Add any labeling desired
    if(! is.null(label_rows)){
        annotation_df$xend1 <- 0.9*annotation_df$xend
        annotation_df$yend1 <- 0.9*annotation_df$yend
        annotation_df$xend2 <- 0.95*annotation_df$xend
        annotation_df$yend2 <- 0.95*annotation_df$yend
        annotation_df$xadj <- 0
        annotation_df$xadj[annotation_df$xend1 < 0] <- 1
        annotation_df$yadj <- 0
        annotation_df$yadj[annotation_df$yend1 < 0] <- 1
        if(colour_code_labels) ac <- annotation_df$cg else ac <- label_colour 
        p <- p + geom_segment(data = annotation_df,
                              aes_string(x = "x",
                                         y = "y",
                                         xend = "xend1",
                                         yend = "yend1"),
                              colour = ac, size = 0.5,
                              arrow = arrow(length = unit(0, "cm"))) +
            geom_text(data = annotation_df,
                      aes_string(x = "xend2", y = "yend2", label = "label"),
                      hjust = "outward",
                      vjust = "outward",
                      color = ac,
                      size = label_size) +
            geom_point(data = annotation_df,
                       aes_string(x = "x", y = "y"),
                       shape = 1,
                       color = "black",
                       size = marker_size)

    }
    return(p)
}
