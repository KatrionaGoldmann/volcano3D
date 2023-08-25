
#' Three-Dimensional Volcano Plot
#'
#' Plots the three-way comparisons of variables such as gene expression data in
#' 3D space using plotly. x, y position represents polar position on 3 axes
#' representing the amount each variable or gene tends to each of the 3
#' categories. The z axis represents -log10 P value for the one-way test
#' comparing each variable across the 3 groups.
#' 
#' @param polar Object of S4 class 'volc3d' following call to either
#'   \code{\link{polar_coords}}, \code{\link{deseq_polar}} or
#'   \code{\link{voom_polar}}
#' @param type Either `1` or `2` specifying type of polar coordinates: `1` =
#'   Z-scaled, `2` = unscaled (equivalent to log2 fold change for gene
#'   expression).
#' @param label_rows A vector of row names or numbers to label
#' @param arrow_length The length of label arrows (default 100)
#' @param label_size font size for labels (default 14).
#' @param colour_code_labels Logical whether label annotations should be colour
#' coded. If `FALSE` `label_colour` is used.
#' @param label_colour HTML colour of annotation labels if not colour coded.
#' @param grid_colour The colour of the cylindrical grid (default "grey80")
#' @param grid_width The width of the grid lines (default 2)
#' @param grid_options Optional list of additional arguments to pass to
#'   \code{\link{polar_grid}}, eg. `z_axis_ticks` and `r_axis_ticks`
#' @param axis_colour The colour of the grid axes and labels (default "black")
#' @param axis_width The width of axis lines (default 2)
#' @param marker_size Size of the markers (default 3)
#' @param marker_outline_width Width for marker outline (default 0 means no
#'   outline)
#' @param marker_outline_colour Colour for marker outline (default white)

#' @param z_axis_title_offset The position scaling between grid and z axis title
#' (default=1.2)
#' @param z_axis_title_size The font size for the z axis title (default=12)
#' @param z_axis_angle Angle in radians for the position of z axis (default
#'   0.5)
#' @param radial_axis_title_size The font size for the radial (default=15)
#' @param radial_axis_title_offset The position scaling between grid and radial 
#' axis title (default=1.2)
#' @param xy_aspectratio The aspect ratio for the xy axis compared to z (default
#'   1). Increasing this makes the grid wider in the plot window.
#' @param z_aspectratio The aspect ratio for the z axis compared to x and y
#'   (default 0.8). Decreasing this makes the plot appear more squat.
#' @param camera_eye The (x,y,z) components of the start 'eye' camera vector.
#'   This vector determines the view point about the origin of this scene.
#' @param ... Optional arguments passed to \code{\link[plotly:plot_ly]{plot_ly}}
#' @return Returns a cylindrical 3D plotly plot featuring variables on a 
#' tri-axis radial graph with the -log10(multi-group test p-value) on the 
#' z-axis
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://pubmed.ncbi.nlm.nih.gov/31461658/}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and 
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @seealso \code{\link{polar_coords}}
#' 
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
#'                           data = t(syn_example_rld))
#' volcano3D(syn_polar)
#' 
#' @importFrom plotly plot_ly add_trace add_text layout
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @export
#' 
volcano3D <- function(polar, type = 1,
                      label_rows = c(),
                      label_size = 14,
                      arrow_length = 100, 
                      colour_code_labels = FALSE,
                      label_colour = "black",
                      grid_colour = "grey80",
                      grid_width = 2,
                      grid_options = NULL,
                      axis_colour = "black",
                      axis_width = 2,
                      marker_size = 3,
                      marker_outline_width = 0,
                      marker_outline_colour = "white",
                      z_axis_title_offset = 1.2,
                      z_axis_title_size = 12,
                      z_axis_angle = 0.5,
                      radial_axis_title_size = 14, 
                      radial_axis_title_offset = 1.2,
                      xy_aspectratio = 1,
                      z_aspectratio = 0.8,
                      camera_eye = list(x=0.9, y=0.9, z=0.9),
                      ...) {
  
  # Check the input data
  if (is(polar, "polar")) {
    args <- as.list(match.call())[-1]
    return(do.call(volcano3D_v1, args))  # for back compatibility
  }
  if (!is(polar, "volc3d")) stop("Not a 'volc3d' class object")
  args <- list(r_vector = polar@df[[type]]$r, z_vector = polar@df[[type]]$z)
  args <- append(args, grid_options)
  
  # Build the coordinate system
  grid <- do.call(polar_grid, args)
  polar_grid <- grid@polar_grid
  axes <- grid@axes
  axis_labels <- grid@axis_labs
  h <- grid@z
  R <- grid@r
  max_offset <- max(c(z_axis_title_offset, radial_axis_title_offset))
  xyrange <- c(-1.05*(max_offset)*R, 
               1.05*max_offset*R)
  axis_settings <- list(title = "", zeroline = FALSE, showline = FALSE, 
                        showticklabels = FALSE, showgrid = FALSE, 
                        autotick = FALSE, showspikes = FALSE)
  axis_settings_xy <- list(title = "", zeroline = FALSE, showline = FALSE, 
                           showticklabels = FALSE, showgrid = FALSE, 
                           autotick = FALSE, showspikes = FALSE,
                           range = xyrange)
  df <- polar@df[[type]]
  
  # Label rows
  if(length(label_rows) != 0){
    if(! all(is.numeric(label_rows))) {
      if(! all(label_rows %in% rownames(df))) {
        stop("label_rows must be in rownames(polar_df)")
      }}
    if(all(is.numeric(label_rows))) {
      if(! all(label_rows < nrow(df))) {
        stop("label_rows not in 1:nrow(polar_df)")
      }}
    annot <- lapply(label_rows, function(i) {
      row <- df[i, ]
      if(colour_code_labels) ac <- row$col else ac <- label_colour
      ac <- col2hex(ac)
      annot <- list(x = row$x, y = row$y, z = row$z, 
                    text = rownames(row), 
                    textangle = 0, ax = arrow_length, ay  = 0,
                    arrowcolor = ac, font = list(color = ac),
                    arrowwidth = 1, arrowhead = 0, arrowsize = 1.5, 
                    yanchor = "middle")
    })
  } else {annot <- list()}
  
  # Generate plotly plot
  plot_ly(polar@df[[type]], x = ~x, y = ~y, z = ~z, color = ~lab, 
          colors = polar@scheme,  hoverinfo='text',
          text = ~paste0(rownames(polar@df[[type]]), 
                         "<br>theta = ", as.integer(angle),
                         ", r = ", formatC(r, digits = 3),
                         "<br>P = ", 
                         format(pvalue, digits = 3, scientific = 3)),
          key = rownames(polar@df[[type]]),
          marker = list(size = marker_size,
                        line = list(color = marker_outline_colour, 
                                    width = marker_outline_width)),
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
      z = 0, text = levels(polar@outcome),
      color = I(axis_colour), type = "scatter3d", mode = "text", 
      textfont = list(size = radial_axis_title_size), 
      textposition = 'middle center', 
      hoverinfo = 'none', showlegend = FALSE, inherit = FALSE) %>%
    
    # label z axis
    add_text(x = c(rep(R*sinpi(z_axis_angle), 
                       grid@n_z_breaks), 
                   R*z_axis_title_offset*sinpi(z_axis_angle)),
             y = c(rep(R*cospi(z_axis_angle), 
                       grid@n_z_breaks), 
                   R*z_axis_title_offset*cospi(z_axis_angle)),
             z = c(grid@z_breaks, h/2),
             text = c(grid@z_breaks, '-log<sub>10</sub>P'),
             textposition = 'middle left', 
             textfont = list(size = z_axis_title_size),  
             color = I(axis_colour), hoverinfo = 'none', 
             showlegend = FALSE, inherit = FALSE) %>%
    
    # label radial axis
    add_text(x = grid@text_coords$x, 
             y = grid@text_coords$y, z = h * 0.03,
             text = grid@text_coords$text, textposition = 'middle center', 
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
}


# from gplots
col2hex <- function(cname) {
  colMat <- grDevices::col2rgb(cname)
  grDevices::rgb(red = colMat[1,] / 255,
                 green = colMat[2,] / 255,
                 blue = colMat[3,] / 255)
}
