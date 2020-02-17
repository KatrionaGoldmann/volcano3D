#' Polar Grid structure for Three Way Polar Plot
#'
#' This function creates polar grid of radius r
#' @param r_vector numerical vector of radial coordinates
#' @param axis_ticks a numerical vector of breaks for the radial axis
#' @param axis_angle Angle to position the radial axis in pi (default = 4/3)
#' @param ... optional parameters for \link[base]{pretty}
#' @return Returns a list containing:
#' \itemize{
#'   \item{'polar_grid'} The coordinates for a radial grid
#'   \item{'axes'} The axes features for plotly
#'   \item{'axis_labels'} The axis labels
#'   \item{'r'} The grid radius
#'   \item{'text_coords'} The coordinates for text labels
#' }
#' @keywords pvalue, polar, plot
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' "Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @export
#' @examples
#' polar_grid(r_vector = 1:10, axis_ticks = 1:10)

polar_grid <- function(r_vector, axis_ticks, axis_angle = 4/3, ...){
  
  if(is.null(axis_ticks)) {
    r_breaks <- pretty(r_vector, ...)
  } else{ r_breaks <- axis_ticks}
  r_breaks <- r_breaks[! is.na(r_breaks)]
  n_breaks <- length(r_breaks) - 1
  
  ## To make the rings continuous we need frequnt dots
  ## cos pi works in pi so we want a dot every 1/100 or a pi
  circular_grid <- data.frame(x = unlist(lapply(1:n_breaks, function(b) {
    c(b*max(r_breaks)/n_breaks*cospi(0:200/100), NA)
  })),
  y = unlist(lapply(1:n_breaks, function(b) {
    c(b*max(r_breaks)/n_breaks*sinpi(0:200/100), NA)}
  )))
  
  # radial spokes out (all grey spokes - not the black ones, 9 total)
  polar_grid_top <- data.frame(x = unlist(lapply(c(1:3, 5:7, 9:11), function(i){
    c(max(r_breaks)/n_breaks*cospi(i*2/12),
      rep(max(r_breaks)*cospi(i*2/12), 2), NA)
  })),
  y = unlist(lapply(c(1:3, 5:7, 9:11), function(i){
    c(max(r_breaks)/n_breaks*sinpi(i*2/12),
      rep(max(r_breaks)*sinpi(i*2/12), 2), NA)
  })))
  polar_grid <- rbind(circular_grid, polar_grid_top)
  
  # The LMF axes
  axes <- data.frame(x = unlist(lapply(0:2, function(i){
    c((max(r_breaks)/n_breaks)*cospi(i*2/3),
      rep(max(r_breaks)*cospi(i*2/3), 2), NA)
  })),
  y = unlist(lapply(0:2, function(i){
    c((max(r_breaks)/n_breaks)*sinpi(i*2/3),
      rep(max(r_breaks)*sinpi(i*2/3), 2), NA)
  })))
  radial_spokes <- data.frame(x = rep(0,3), y = rep(0,3),
                              xend = cospi(0:2 * 2/3), yend = sinpi(0:2 * 2/3))
  
  axis_labs <- data.frame(x = radial_spokes$xend*max(r_breaks),
                          y = radial_spokes$yend*(max(r_breaks)))
  
  axis_labs$x_adjust <- unlist(lapply(sign(axis_labs$x), function(s) {
    switch(as.character(s), "1" = "right", "-1" = "left", "0" = "center")
  }))
  axis_labs$y_adjust <- unlist(lapply(sign(axis_labs$y), function(s) {
    switch(as.character(s), "1" = "top", "-1" = "bottom", "0" = "middle")
  }))
  axis_labs$adjust <- paste(axis_labs$y_adjust, axis_labs$x_adjust)
  
  text_coords <- data.frame(x = r_breaks[2:length(r_breaks)]*sinpi(axis_angle),
                           y = r_breaks[2:length(r_breaks)]*cospi(axis_angle),
                           text = format(r_breaks[2:length(r_breaks)], 
                                         digits = 2))
  
  return(list("polar_grid" = polar_grid,
              "axes" = 	axes,
              "axis_labs" = axis_labs,
              "r" = max(r_breaks),
              "text_coords" = text_coords))
}

#' Grid required for 3D volcano plot
#'
#' Generates a cylindrical grid of the appropriate dimensions for a 3D volcano
#' plot
#' @param r_vector numerical vector of radial coordinates
#' @param r_axis_ticks a numerical vector of breaks for the radial axis
#' @param z_vector numerical vector of z coordinates
#' @param z_axis_ticks a numerical vector of breaks for the z axis
#' @param axis_angle angle to position the radial axis in pi (default = 5/6)
#' @param ... optional parameters for \code{\link[base]{pretty}}
#' @return Returns a list containing:
#' \itemize{
#'   \item{'polar_grid'} The coordinates for a radial grid
#'   \item{'axes'} The axes features for plotly
#'   \item{'axis_labels'} The axis labels
#'   \item{'r'} The grid radius
#'   \item{'z'} The grid height
#'   \item{'text_coords'} The coordinates for text labels
#'   \item{'n_r_breaks'} The number of ticks on the r axis
#'   \item{'n_r_breaks'} The number of ticks on the z axis
#' }
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @keywords pvalue, polar, plot
#' @export
#' @examples
#' polar_grid3D(r_vector = 1:10, 
#'            z_vector = 1:5, 
#'            r_axis_ticks = 1:10, 
#'            z_axis_ticks = 1:5)

polar_grid3D <- function(r_vector,
                         z_vector,
                         r_axis_ticks,
                         z_axis_ticks,
                         axis_angle = 5/6, 
                         ...){
  
  if(is.null(r_axis_ticks)) {
    r_breaks <- pretty(r_vector, ...)
  } else{ r_breaks <- r_axis_ticks}
  r_breaks <- r_breaks[! is.na(r_breaks)]
  
  if(is.null(z_axis_ticks)) {
    z_breaks <- pretty(z_vector)
  } else{ z_breaks <- z_axis_ticks}
  z_breaks <- z_breaks[! is.na(z_breaks)]
  
  n_r_breaks <- length(r_breaks) - 1
  n_z_breaks <- length(z_breaks) - 1
  
  # Set up the concentric circles on the x/y plane 
  # (Circles split by NA to make discontinuous)
  cylindrical_grid <- data.frame(x = unlist(lapply(1:n_r_breaks, function(i){
    c(max(r_breaks)/n_r_breaks*i*cospi(0:200/100), NA)
  })),
  y = unlist(lapply(1:n_r_breaks, function(i){
    c(max(r_breaks)/n_r_breaks*i*sinpi(0:200/100), NA)
  })),
  z = 0)
  
  # radial spokes out
  polar_grid_top <- data.frame(x = unlist(lapply(c(1:3, 5:7, 9:11), function(i){
    c(max(r_breaks)/n_r_breaks*cospi(i*2/12),
      rep(max(r_breaks)*cospi(i*2/12), 2), NA)
  })),
  y = unlist(lapply(c(1:3, 5:7, 9:11), function(i){
    c(max(r_breaks)/n_r_breaks*sinpi(i*2/12),
      rep(max(r_breaks)*sinpi(i*2/12), 2), NA)
  })),
  z = rep(c(0,0,max(z_breaks),NA), 9))
  
  # Create the circles on the cylinder - h/d cylinders
  z_cyl <- c()
  for (i in (1:n_z_breaks)*(max(z_breaks)/n_z_breaks)){
    z_cyl <- c(z_cyl, rep(i, 201), NA)
  }
  cylinder <- data.frame(x = rep(c(max(r_breaks)*cospi(0:200/100), NA), 
                                 n_z_breaks),
                         y = rep(c(max(r_breaks)*sinpi(0:200/100), NA), 
                                 n_z_breaks),
                         z = z_cyl)
  
  polar_grid <- rbind(polar_grid_top, cylindrical_grid, cylinder)
  
  # Add the three axes
  axes <- data.frame(x = unlist(lapply(0:2, function(i){
    c(max(r_breaks)/n_r_breaks*cospi(i*2/3), 
      rep(max(r_breaks)*cospi(i*2/3), 2), NA)
  })),
  y = unlist(lapply(0:2, function(i){
    c(max(r_breaks)/n_r_breaks*sinpi(i*2/3), 
      rep(max(r_breaks)*sinpi(i*2/3), 2), NA)
  })),
  z = rep(c(0, 0, max(z_breaks), NA), 3))
  radial_spokes <- data.frame(x = rep(0,3),
                              y = rep(0,3),
                              xend = cospi(0:2 * 2/3),
                              yend = sinpi(0:2 * 2/3))
  
  axis_labs <- data.frame(x = radial_spokes$xend*max(r_breaks),
                          y = radial_spokes$yend*(max(r_breaks)) )
  
  axis_labs$x_adjust <- unlist(lapply(sign(axis_labs$x), function(s) {
    switch(as.character(s), "1" = "right", "-1" = "left", "0" = "center")
  }))
  axis_labs$y_adjust <- unlist(lapply(sign(axis_labs$y), function(s) {
    switch(as.character(s), "1" = "top", "-1" = "bottom", "0" = "middle")
  }))
  axis_labs$adjust <- paste(axis_labs$y_adjust, axis_labs$x_adjust)
  
  text_coords <- data.frame(x = r_breaks[2:length(r_breaks)]*sinpi(axis_angle),
                           y = r_breaks[2:length(r_breaks)]*cospi(axis_angle),
                           text = format(r_breaks[2:length(r_breaks)], 
                                         digits = 2))
  
  
  return(list("polar_grid" = polar_grid,
              "axes" = axes,
              "axis_labs" = axis_labs,
              "r" = max(r_breaks),
              "z" = max(z_breaks),
              "text_coords" = text_coords,
              "n_r_breaks" = n_r_breaks,
              "n_z_breaks" = n_z_breaks))
}
