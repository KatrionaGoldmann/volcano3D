

#' Three-Dimensional Volcano Plot
#'
#' Version 2.0. Plots the three-way comparisons of variables such as gene
#' expression data in 3D space using plotly. x, y position represents polar
#' position on 3 axes representing the amount each variable or gene tends to
#' each of the 3 categories. The z axis represents -log10 P value for the
#' one-way test comparing each variable across the 3 groups.
#' 
#' @param obj Object of S3 class 'volc3d' following call to either
#'   `polarCoords()` or `DESeqToVolc()`
#' @param type Either `1` or `2` specifying type of polar coordinates: `1` =
#'   Z-scaled, `2` = unscaled (equivalent to log2 fold change for gene
#'   expression).
#' @param grid_colour The colour of the cylindrical grid (default="grey80"). 
#' @param axis_colour The colour of the grid axes and labels (default="black").
#' @param marker_size Size of the markers (default = 2.7).
#' @param marker_outline_colour Colour for marker outline (default = white)
#' @param marker_outline_width Width for marker outline (default = 0)
#' @param axis_angle Angle in radians for the position of z axis (default 0.5). 
#' @param z_aspectratio The aspect ratio for the z axis compared to x and y.
#'   Decreasing this makes the plot appear more squat. If `NULL` it is set
#'   automatically.
#' @param xy_aspectratio The aspect ratio for the xy axis compared to z (default
#'   = 1). Increasing this makes the grid wider in the plot window.
#' @param camera_eye The (x,y,z) components of the start 'eye' camera vector.
#'   This vector determines the view point about the origin of this scene.
#' @importFrom plotly plot_ly add_trace add_text layout %>%
#' @export
#' 
volcano3dx <- function(obj, type = 1,
                       axis_width = 2,
                       grid_colour = "grey80",
                       axis_colour = "black",
                       grid_width = 2,
                       z_axis_title_offset = 1.2,
                       z_axis_title_size = 12,
                       radial_axis_title_size = 14, 
                       radial_axis_title_offset = 1.2,
                       axis_angle = 0.5,
                       xy_aspectratio = 0.75,
                       z_aspectratio = NULL,
                       camera_eye = list(x=0.75, y=0.75, z=0.75),
                       ...) {
  if (!inherits(obj, "volc3d")) stop("Not a 'volc3d' class object")
  grid <- polar_grid(obj[[type]]$r, obj[[type]]$z)
  polar_grid <- grid@polar_grid
  axes <- grid@axes
  axis_labels <- grid@axis_labs
  h <- grid@z
  R <- grid@r
  if (is.null(z_aspectratio)) z_aspectratio <- 18 / h
  xyrange <- c(-1.05*z_axis_title_offset*R, 
               1.05*z_axis_title_offset*R)
  axis_settings <- list(title = "", zeroline = FALSE, showline = FALSE, 
                        showticklabels = FALSE, showgrid = FALSE, 
                        autotick = FALSE, spikesides = FALSE)
  axis_settings_xy <- list(title = "", zeroline = FALSE, showline = FALSE, 
                           showticklabels = FALSE, showgrid = FALSE, 
                           autotick = FALSE, showspikes = FALSE,
                           range = xyrange)
  
  plot_ly(obj[[type]], x = ~x, y = ~y, z = ~z, color = ~lab, colors = obj$scheme,
          hoverinfo='text',
          text = ~paste0(rownames(obj[[type]]), "<br>theta = ", as.integer(angle),
                         "<br>r = ", formatC(r, digits = 3),
                         "<br>P = ", format(pvalue, digits = 3, scientific = 3)),
          marker = list(size = 3),
          type = "scatter3d", mode = "markers", ...) %>%
    
    # Add the cylindrical grid
    add_trace(x = polar_grid$x, y = polar_grid$y, z = polar_grid$z, 
              color = I(grid_colour), line = list(width = grid_width),
              showlegend = FALSE, type = "scatter3d", mode = "lines", 
              hoverinfo = "none", inherit = FALSE) %>%
    
    # Axes
    add_trace(x = axes$x, y = axes$y, z = axes$z, color = I(axis_colour),
              line = list(width = axis_width), showlegend = FALSE, 
              type = "scatter3d", mode = "lines", hoverinfo = "none", 
              inherit = FALSE) %>%
    
    add_text(
      x = radial_axis_title_offset*axis_labels$x, 
      y = radial_axis_title_offset*axis_labels$y, 
      z = 0, text = levels(obj$outcome),
      color = I(axis_colour), type = "scatter3d", mode = "text", 
      textfont = list(size = radial_axis_title_size), 
      textposition = 'middle center', 
      hoverinfo = 'none', showlegend = FALSE, inherit = FALSE) %>%
    
    # label z axis
    add_text(x = c(rep(R*sinpi(axis_angle), 
                       grid@n_z_breaks), 
                   R*z_axis_title_offset*sinpi(axis_angle)),
             y = c(rep(R*cospi(axis_angle), 
                       grid@n_z_breaks), 
                   R*z_axis_title_offset*cospi(axis_angle)),
             z = c(grid@z_breaks, h/2),
             text = c(grid@z_breaks, '-log<sub>10</sub>P'),
             textposition = 'middle left', 
             textfont = list(size = z_axis_title_size),  
             color = I(axis_colour), hoverinfo = 'none', 
             showlegend = FALSE, inherit = FALSE) %>%
    
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
        zaxis = axis_settings
      ),
      xaxis = list(title = "x"),
      yaxis = list(title = "y")
    )
}

