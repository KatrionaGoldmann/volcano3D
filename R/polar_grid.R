setClassUnion("numeric_or_integer_or_NULL", c("numeric", "integer", "NULL"))

#' An S4 class to define the polar grid coordinates system.
#' 
#' @slot polar_grid The coordinates for the cylindrical grid segments with 
#' x,y,z coordinates
#' @slot axes The axes features for plotly
#' @slot axis_labs The axis labels
#' @slot r The grid radius
#' @slot z The grid height
#' @slot text_coords data frame for axis label cartesian coordinates (x, y, z)
#' @slot n_r_breaks The number of ticks on the r axis
#' @slot n_z_breaks The number of ticks on the z axis
#' @slot r_breaks The r axis ticks as a numeric
#' @slot z_breaks The z axis ticks as a numeric

setClass("grid", slots = list(
  polar_grid = "data.frame",
  axes = "data.frame",
  axis_labs = "list",
  r = "numeric",
  z = "numeric",
  text_coords = "data.frame",
  n_r_breaks = "numeric_or_integer_or_NULL",
  n_z_breaks = "numeric_or_integer_or_NULL", 
  r_breaks = "numeric_or_integer_or_NULL",
  z_breaks = "numeric_or_integer_or_NULL"
))

#' Grid required for 3D volcano plot and 2D radial plots
#'
#' Generates a cylindrical grid of the appropriate dimensions for a 3D volcano
#' plot
#' @param r_vector An optional numerical vector for the  radial coordinates. 
#' This is used to calculate breaks on the r axis using 
#' \code{\link[base]{pretty}}. If this is NULL the r_axis_ticks are used as 
#' breaks. 
#' @param z_vector An optional numerical vector for the z coordinates. 
#' This is used to calculate breaks on the z axis using \code{pretty}. If this 
#' is NULL the z_axis_ticks are used as breaks. 
#' @param r_axis_ticks A numerical vector of breaks for the radial axis (used 
#' if r_vector is NULL). 
#' @param z_axis_ticks A numerical vector of breaks for the z axis (used 
#' if z_vector is NULL). 
#' @param axis_angle angle to position the radial axis in pi (default = 5/6)
#' @param n_spokes the number of outward spokes to be plotted (default = 12)
#' @param ... optional parameters for \code{\link[base]{pretty}} on the r axis
#' @return Returns an S4 grid object containing:
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
#' library(volcano3Ddata)
#' data(syn_data)
#' syn_p_obj <- create_dep(sampledata = syn_metadata, 
#'                     contrast = "Pathotype", 
#'                     pvalues = syn_pvalues,
#'                     p_col_suffix = "pvalue", 
#'                     fc_col_suffix = "log2FoldChange",
#'                     multi_group_prefix = "LRT", 
#'                     expression = syn_rld)
#'                     
#' polar_grid(r_vector=syn_polar@polar$r_zscore,
#'            z_vector=NULL,
#'            r_axis_ticks = NULL,
#'            z_axis_ticks = c(0, 8, 16, 32),
#'            n_spokes = 4)

polar_grid <- function(r_vector = NULL,
                       z_vector = NULL,
                       r_axis_ticks = NULL,
                       z_axis_ticks = NULL,
                       axis_angle = 5/6, 
                       n_spokes = 12, 
                       ...){
  
  if(is.null(r_axis_ticks)) {
    r_breaks <- pretty(r_vector, ...)
  } else{ r_breaks <- r_axis_ticks}
  r_breaks <- r_breaks[! is.na(r_breaks)]
  r_breaks <- sort(r_breaks)
  if(r_breaks[1] != 0) r_breaks <- c(0, r_breaks)
  
  if(is.null(z_axis_ticks)) {
    z_breaks <- pretty(z_vector)
  } else{ z_breaks <- z_axis_ticks }
  z_breaks <- z_breaks[! is.na(z_breaks)]
  z_breaks <- sort(z_breaks)
  if(length(z_breaks) > 0) { if( z_breaks[1] != 0) z_breaks <- c(0, z_breaks)}
  
  n_r_breaks <- length(r_breaks) - 1
  n_z_breaks <- length(z_breaks) - 1
  if(n_z_breaks < 0) { n_z_breaks <- 0 }
  
  # Set up the concentric circles on the x/y plane 
  # (Circles split by NA to make discontinuous)
  cylindrical_grid <- data.frame(
    x = unlist(lapply(1:n_r_breaks, function(i){
      c(max(r_breaks)/n_r_breaks*i*cospi(0:200/100), NA)
    })),
    y = unlist(lapply(1:n_r_breaks, function(i){
      c(max(r_breaks)/n_r_breaks*i*sinpi(0:200/100), NA)
    })),
    z = 0, area = "cylindrical_grid")
  
  # radial spokes out
  mz <- switch(as.character(is.null(z_breaks)), "TRUE"=0, "FALSE"=max(z_breaks))
  
  polar_grid_top <- data.frame(
    x = unlist(lapply(c(1:n_spokes), function(i){
      c(max(r_breaks)/n_r_breaks*cospi(i*2/n_spokes),
        rep(max(r_breaks)*cospi(i*2/n_spokes), 2), NA)
    })),
    y = unlist(lapply(c(1:n_spokes), function(i){
      c(max(r_breaks)/n_r_breaks*sinpi(i*2/n_spokes),
        rep(max(r_breaks)*sinpi(i*2/n_spokes), 2), NA)
    })),
    z = rep(c(0, 0, mz, NA), n_spokes), 
    area = "polar grid top")
  
  # Create the circles on the cylinder - h/d cylinders
  z_cyl <- c()
  for (i in z_breaks[2:length(z_breaks)]){
    z_cyl <- c(z_cyl, rep(i, 201), NA)
  }
  
  if(is.null(z_vector) & is.null(z_axis_ticks)){
    cylinder <- NULL
  } else{
    cylinder <- data.frame(x = rep(c(max(r_breaks)*cospi(0:200/100), NA), 
                                   n_z_breaks),
                           y = rep(c(max(r_breaks)*sinpi(0:200/100), NA), 
                                   n_z_breaks),
                           z = z_cyl, 
                           area = "cylinder")
  }
  polar_grid <- rbind(polar_grid_top, cylindrical_grid, cylinder)
  
  # Add the three axes
  axes <- data.frame(
    x = unlist(lapply(0:2, function(i){
      c(max(r_breaks)/n_r_breaks*cospi(i*2/3), 
        rep(max(r_breaks)*cospi(i*2/3), 2), NA)
    })),
    y = unlist(lapply(0:2, function(i){
      c(max(r_breaks)/n_r_breaks*sinpi(i*2/3), 
        rep(max(r_breaks)*sinpi(i*2/3), 2), NA)
    })),
    z = rep(c(0, 0, mz, NA), 3))
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
  
  
  methods::new("grid",
              polar_grid = polar_grid,
              axes = axes,
              axis_labs = axis_labs,
              r = max(r_breaks),
              z = mz,
              text_coords = text_coords,
              n_r_breaks = n_r_breaks,
              n_z_breaks = n_z_breaks, 
              r_breaks = r_breaks[2:length(r_breaks)],
              z_breaks = z_breaks[2:length(z_breaks)])
}
