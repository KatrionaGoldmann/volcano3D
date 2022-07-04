
#' @importFrom grDevices col2rgb

radial_ggplot_v1 <- function(polar,
                          colours = c("green3", "cyan", "blue",
                                      "purple", "red", "gold2"),
                          non_sig_colour = "grey60",
                          colour_scale = "discrete",
                          continuous_shift = 1.33,
                          label_rows = NULL,
                          arrow_length = 1,
                          grid = NULL,
                          fc_or_zscore = "zscore",
                          label_size = 5,
                          colour_code_labels = TRUE,
                          label_colour = NULL,
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
                          axis_angle = 1/6,
                          legend_size = 20,
                          ...){

    if(! class(polar) %in% c("polar")) stop("polar must be a polar object")
    polar_df <- polar@polar

    if(! is.null(colours) & length(colours) != 6){
        stop(paste("colours must be a character vector of plotting colours of",
                   "length six. One colour for each significance group"))
    }
    if(! is.numeric(continuous_shift)) {
        stop('continuous_shift must be numeric')
    }
    if(! (0 <= continuous_shift & continuous_shift <= 2) ) {
        stop('continuous_shift must be between 0 and 2')
    }
    if(! colour_code_labels & is.null(label_colour)){
        stop('If colour_code_labels is false please enter a valid label_colour')
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
    colours <- unlist(lapply(c(colours, non_sig_colour), function(x) {
        if(! grepl("#", x) &
           class(try(col2rgb(x), silent = TRUE))[1] == "try-error") {
            stop(paste(x, 'is not a valid colour'))
        } else if (! grepl("#", x) ) {
            y <- col2rgb(x)[, 1]
            x <- rgb(y[1], y[2], y[3], maxColorValue=255)
        }
        return(x)
    }))
    colours <- setNames(colours, c(sig_groups, polar@non_sig_name))

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

    # Calculate the continuous colours
    offset <- (polar_df$angle[!is.na(polar_df$angle)] + continuous_shift/2)
    offset[offset > 1] <- offset[offset > 1] - 1
    polar_df$hue <- hsv(offset, 1, 1)
    polar_df$hue[polar_df$sig == polar@non_sig_name] <- non_sig_colour

    # account for duplicated colours
    if(any(duplicated(colours))){
        warning(paste("Some colours are repeated. These will be compressed",
                      "into one significance group"))

        colours <- setNames(unique(colours),
                            unlist(lapply(unique(colours), function(x) {
                                paste(names(colours)[colours == x],
                                      collapse=" or \n")
                            })))
        polar_df$sig <- factor(names(colours)[match(polar_df$col, colours)])
    }

    # make sure the non-sig markers are on the bottom - reshuffle the order
    polar_df$sig <-
        factor(polar_df$sig,
               levels = c(polar@non_sig_name,
                          as.character(
                              unique(polar_df$sig[polar_df$sig !=
                                                      polar@non_sig_name]))))

    cols <- colours[match(levels(droplevels(polar_df$sig)), names(colours))]

    if(is.null(grid)) {
        grid <- polar_grid(r_vector = polar_df$r,
                           r_axis_ticks = NULL,
                           axis_angle = axis_angle,
                           ...)
    } else {  if(class(grid) != "grid") stop('grid must be a grid object')}

    grid@polar_grid <- grid@polar_grid[grid@polar_grid$area != "cylinder", ]

    # markers are plotted in order of rows so push ns to the bottom of plot
    polar_df <- polar_df[c(which(polar_df$sig ==  polar@non_sig_name),
                           which(polar_df$sig !=  polar@non_sig_name)), ]

    # alignment for text
    hadj <- -1*sign(grid@axis_labs$x)
    hadj[hadj ==  -1] <- 0

    polar_df$cg <- switch(colour_scale,
                          "discrete"=polar_df$col,
                          "continuous"=polar_df$hue)

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
                 label = levels(polar@sampledata[, polar@contrast]),
                 color = axis_colour, size = axis_title_size) +

        # Add markers
        geom_point(aes_string(fill=switch(colour_scale,
                                            "discrete"= "sig",
                                            "continuous" = "hue")),
                   size = marker_size,
                   alpha = marker_alpha, 
                   color = marker_outline_colour,
                   stroke = marker_outline_width,
                   shape = 21) +

        scale_fill_manual(name="",
            values = switch(colour_scale,
                       "discrete" = as.character(cols),
                       "continuous" = levels(factor(polar_df$hue)))) +

        # Set the background colour etc.
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.key = element_blank(),
              legend.justification = c(1, 1),
              legend.position = switch(colour_scale,
                                       "discrete"="right",
                                       "continuous"="none"),
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
