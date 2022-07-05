setClassUnion("character_or_NULL", c("character", "NULL"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))

setClass("polar", slots = list(sampledata = "data.frame",
                               contrast = "character",
                               pvalues = "data.frame",
                               multi_group_test = "character_or_NULL",
                               expression = "df_or_matrix",
                               polar = "df_or_matrix",
                               non_sig_name = "character"))


volcano3D_v1 <- function(polar,
                      colours=c("green3", "cyan", "blue", 
                                "purple", "red", "gold2"), 
                      non_sig_colour = "grey60",
                      colour_scale = "discrete",
                      continuous_shift = 1.33, 
                      z_axis_title_offset = 1.2,
                      radial_axis_title_offset = 1,
                      z_axis_title_size = 15,
                      radial_axis_title_size = 15, 
                      continue_radial_axes = TRUE,
                      axis_width = 2,
                      grid_width = 2,
                      label_rows = c(),
                      grid = NULL, 
                      fc_or_zscore = "zscore",
                      label_size = 14,
                      arrow_length=50, 
                      colour_code_labels = TRUE,
                      label_colour = NULL,
                      hover_text = "label",
                      grid_colour = "grey80", 
                      axis_colour = "black",
                      marker_size = 2.7,
                      marker_alpha = 1,
                      marker_outline_colour = "white",
                      marker_outline_width = 0,
                      axis_angle = 0.5, 
                      z_aspectratio = 1, 
                      xy_aspectratio = 1,
                      plot_height = 700,
                      camera_eye = list(x=1.25, y=1.25, z=1.25),
                      source="volcano3D",
                      ...){
    
    if(! is(polar, "polar")) stop("polar must be a polar object")
    polar_df <- cbind(polar@polar, polar@pvalues)
    
    polar_df$hover <- eval(parse(
        text = paste0("with(polar_df, ", hover_text, ")")))
    
    polar_df <- polar_df[, c(colnames(polar@polar), "hover")]
    
    if(! colour_code_labels & is.null(label_colour)){
        stop('If colour_code_labels is false please enter a valid label_colour')
    }
    
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
           inherits(try(col2rgb(x), silent = TRUE), "try-error")) {
            stop(paste(x, 'is not a valid colour'))
        } else if (! grepl("#", x) ) {
            y <- col2rgb(x)[, 1]
            x <- rgb(y[1], y[2], y[3], maxColorValue=255)
        }
        return(x)
    }))
    colours <- setNames(colours, c(polar@non_sig_name, sig_groups))
    
    if(! inherits(polar_df, "data.frame")) {
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
    if(! (0 <= continuous_shift & continuous_shift <= 2) ) {
        stop('continuous_shift must be between 0 and 2')
    }
    
    # Calculate the colours by significance
    offset <- (polar_df$angle[!is.na(polar_df$angle)] + continuous_shift/2)
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
    } else { if(inherits(grid, "grid")) stop('grid must be a grid object')}
    
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
            if(colour_code_labels) ac <- row$col else ac <- label_colour 
            annot <- list(x = row$x, y = row$y, z = row$logP, text = row$label, 
                          textangle = 0, ax = arrow_length, ay  = 0,
                          arrowcolor = ac, font = list(color = ac),
                          arrowwidth = 1, arrowhead = 6, arrowsize = 1.5, 
                          yanchor = "middle")
        })
    } else {annot <- list()}
    
    axis_settings_xy[['range']] <- c(-1.05*z_axis_title_offset*(grid@r+1), 
                                     1.05*z_axis_title_offset*(grid@r+1))
    
    p <- plot_ly(data = volcano_toptable, x = ~x, y = ~y, z = ~logP,
                 marker = list(size = marker_size, sizemode = 'diameter', 
                               line = list(color = marker_outline_colour, 
                                           width = marker_outline_width)),
                 height = plot_height,
                 key=~label,
                 opacity = marker_alpha,
                 color = ~switch(colour_scale,
                                 "discrete" = sig,
                                 "continuous" = I(hue)),
                 hoverinfo = 'text', text = ~hover,
                 colors = switch(colour_scale,
                                 "discrete" = colours,
                                 "continuous" = NULL),
                 type = "scatter3d", mode = "markers", source = source) %>%
        
        # Add the cylindrical grid
        add_trace(x = polar_grid$x, y = polar_grid$y, z = polar_grid$z, 
                  color = I(grid_colour), line = list(width = grid_width),
                  showlegend = FALSE, type = "scatter3d", mode = "lines", 
                  hoverinfo = "none",inherit = FALSE) %>%
        
        # Horizontal axes
        add_trace(x = axes$x, y = axes$y, z = 0, color = I(axis_colour),
                  line = list(width = axis_width), showlegend = FALSE, 
                  type = "scatter3d", mode = "lines", hoverinfo = "none", 
                  inherit = FALSE) %>%
        add_text(
            x = radial_axis_title_offset*axis_labels$x*z_axis_title_offset, 
            y = radial_axis_title_offset*axis_labels$y*z_axis_title_offset, 
            z = 0, text = levels(polar@sampledata[, polar@contrast]),
            color = I(axis_colour), type = "scatter3d", mode = "text", 
            textfont = list(size = radial_axis_title_size), 
            textposition = 'middle center', 
            hoverinfo = 'none', showlegend = FALSE, inherit = FALSE) %>%
        
        # label z axis
        add_text(x = c(rep(1.05*R*sinpi(axis_angle), 
                           grid@n_z_breaks), 
                       1.2*R*z_axis_title_offset*sinpi(axis_angle)),
                 y = c(rep(1.05*R*cospi(axis_angle), 
                           grid@n_z_breaks), 
                       1.2*R*z_axis_title_offset*cospi(axis_angle)),
                 z = c(grid@z_breaks, h/2)*0.95,
                 text = c(grid@z_breaks, '-log<sub>10</sub>P'),
                 textposition = 'middle left', 
                 textfont = list(size = z_axis_title_size),  
                 color = I(axis_colour), hoverinfo = 'none', 
                 showlegend = FALSE, inherit = FALSE) %>%
        
        add_trace(x = R*sinpi(axis_angle), y = R*cospi(axis_angle), 
                  z = c(0, h), color = I(axis_colour),
                  line = list(width = axis_width), 
                  showlegend = FALSE, type = "scatter3d", mode = "lines",
                  hoverinfo = "none", inherit = FALSE) %>%
        
        # label radial axis
        add_text(x = grid@text_coords$x, 
                 y = grid@text_coords$y, z = 0.05,
                 text = grid@text_coords$text, textposition = 'top center', 
                 textfont = list(size = 10), color = I(axis_colour), 
                 hoverinfo = 'none', showlegend = FALSE, inherit = FALSE) %>%
        
        layout(
            margin = list(0, 0, 0, 0),
            paper_bgcolor = 'rgba(0, 0, 0, 0)',
            plot_bgcolor = 'rgba(0, 0, 0, 0)',
            
            scene = list(
                camera = list(eye = camera_eye),
                aspectratio = list(x = xy_aspectratio,
                                   y = xy_aspectratio,
                                   z = z_aspectratio),
                dragmode = "turntable",
                xaxis = axis_settings_xy,
                yaxis = axis_settings_xy,
                zaxis = axis_settings,
                annotations = annot
            ),
            xaxis = list(title = "x"),
            yaxis = list(title = "y")
        )
    
    if(continue_radial_axes){
        vals <- c(0, 0, NA, 
                  2*pi/3, 2*pi/3, NA, 
                  4*pi/3, 4*pi/3, NA)
        p <-p %>% add_trace(x = R*sin(pi/2 + vals),
                            y = R*cos(pi/2 + vals),
                            z = rep(c(0, h, NA), 3),
                            color = I(axis_colour), 
                            line = list(width = axis_width),
                            showlegend = FALSE, 
                            type = "scatter3d", mode = "lines", 
                            hoverinfo = "none", inherit = FALSE) 
    }
    return(p)
}

