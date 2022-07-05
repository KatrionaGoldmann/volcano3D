#' Plots grid objects for inspection using plotly
#'
#' This function creates an interactive grids in polar and cylindrical
#' coordinates
#' @param grid A grid object produced by \code{\link{polar_grid}}.
#' @param plot_height The plot height in px (default=700),
#' @param axis_angle The angle in radians at which to add axis (default=0). 
#' @param z_axis_title_offset Offset for z axis title (default=1.2).
#' @return Returns a list containing a polar and cylindrical coordinate system.
#' @importFrom plotly plot_ly add_trace add_text layout %>%
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://doi.org/10.1016/j.celrep.2019.07.091}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @keywords hplot iplot
#' @export
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
#'                           data = t(syn_example_rld))
#'
#' grid <- polar_grid(r_vector=syn_polar@df[[1]]$r,
#'         z_vector=syn_polar@df[[1]]$z,
#'         r_axis_ticks = NULL,
#'         z_axis_ticks = NULL)
#' p <- show_grid(grid)
#' p$polar
#' p$cyl

show_grid <- function(grid, 
                      plot_height=700, 
                      axis_angle=0, 
                      z_axis_title_offset=1.2){
    if(! is.numeric(axis_angle)) {
        stop('axis_angle must be a numeric')
    }
    if(! is.numeric(plot_height)) {
        stop('plot_height must be a numeric')
    }
    if(class(grid)[1] != "grid") {
        stop('grid must be a grid object created by polar_grid()')
    }

    axis_settings <- list(title = "", zeroline = FALSE, showline = FALSE,
                          showticklabels = FALSE, showgrid = FALSE,
                          autotick = FALSE)
    axis_settings_xy <- list(title = "", zeroline = FALSE, showline = FALSE,
                             showticklabels = FALSE, showgrid = FALSE,
                             autotick = FALSE, showspikes = FALSE)
    axis_settings_xy[['range']] <- c(-grid@r, grid@r) * 1.05 * z_axis_title_offset

    grid_cyl <- grid
    cyl <- plot_ly(height = plot_height) %>%
        add_trace(x = grid_cyl@polar_grid$x,
                  y = grid_cyl@polar_grid$y,
                  z = grid_cyl@polar_grid$z,
                  color = I("#CBCBCB"),
                  line = list(width = 2),
                  showlegend = FALSE,
                  type = "scatter3d",
                  mode = "lines",
                  hoverinfo = "none",
                  inherit = FALSE) %>%
        # Horizontal axes
        add_trace(x = grid_cyl@axes$x,
                  y = grid_cyl@axes$y,
                  z = 0,
                  color = I("black"),
                  line = list(width = 2),
                  showlegend = FALSE,
                  type = "scatter3d",
                  mode = "lines",
                  hoverinfo = "none",
                  inherit = FALSE) %>%
        # Axis labels
        add_text(x = grid_cyl@axis_labs$x *1.1,
                 y = grid_cyl@axis_labs$y *1.1,
                 z = 0,
                 text = c("A", "B", "C"),
                 color = I("black"),
                 type = "scatter3d",
                 mode = "text",
                 textfont = list(size = 16),
                 textposition = 'middle center',
                 hoverinfo = 'none',
                 showlegend = FALSE,
                 inherit = FALSE) %>%
        # label z axis
        add_text(x = c(rep(1.05*grid_cyl@r*sinpi(axis_angle),
                           grid_cyl@n_z_breaks),
                       z_axis_title_offset*grid_cyl@r*sinpi(axis_angle)),
                 y = c(rep(1.05*grid_cyl@r*cospi(axis_angle),
                           grid_cyl@n_z_breaks),
                       z_axis_title_offset*grid_cyl@r*cospi(axis_angle)),
                 z = c(grid_cyl@z_breaks, grid_cyl@z/2)*0.95,
                 text = c(grid_cyl@z_breaks, '-log<sub>10</sub>P'),
                 textposition = 'middle left',
                 textfont = list(size = 10),
                 color = I("black"),
                 hoverinfo = 'none',
                 showlegend = FALSE,
                 inherit = FALSE) %>%

        add_trace(x = grid_cyl@r*sinpi(axis_angle),
                  y = grid_cyl@r*cospi(axis_angle),
                  z = c(0, grid_cyl@z),
                  color = I("black"),
                  line = list(width = 2),
                  showlegend = FALSE,
                  type = "scatter3d",
                  mode = "lines",
                  hoverinfo = "none",
                  inherit = FALSE) %>%
        # label radial axis
        add_text(x = grid_cyl@text_coords$x,
                 y = grid_cyl@text_coords$y, z = grid_cyl@z * 0.03,
                 text = grid_cyl@text_coords$text,
                 textposition = 'middle center',
                 textfont = list(size = 10),
                 color = I("black"),
                 hoverinfo = 'none',
                 showlegend = FALSE,
                 inherit = FALSE) %>%

        layout(showlegend = TRUE,
               dragmode = "turntable",
               margin = list(0, 0, 0, 0),
               paper_bgcolor = 'rgba(0, 0, 0, 0)',
               plot_bgcolor = 'rgba(0, 0, 0, 0)',
               scene = list(xaxis = axis_settings_xy,
                            yaxis = axis_settings_xy,
                            zaxis = axis_settings,
                            aspectratio = list(x = 1,
                                               y = 1,
                                               z = 0.8)))

    grid_polar <- grid
    grid_polar@polar_grid <- grid_polar@polar_grid[grid_polar@polar_grid$area
                                                   != "cylinder", ]

    polar <- plot_ly(height = plot_height) %>%
        # add the grid
        add_trace(x = grid_polar@polar_grid$x,
                  y = grid_polar@polar_grid$y,
                  color = I("#CBCBCB"),
                  line = list(width = 1),
                  showlegend = FALSE,
                  type = "scatter",
                  mode = "lines",
                  hoverinfo = "none") %>%
        # add the "horizontal" axes
        add_trace(x = grid_polar@axes$x,
                  y = grid_polar@axes$y,
                  color = I("black"),
                  line = list(width = 2),
                  showlegend = FALSE,
                  type = "scatter",
                  mode = "lines",
                  hoverinfo = "none",
                  inherit = FALSE) %>%
        # add the label text
        add_text(x = grid_polar@axis_labs$x,
                 y = grid_polar@axis_labs$y,
                 text = c("A", "B", "C"),
                 color = I("black"),
                 type = "scatter",
                 mode = "text",
                 textfont = list(size = 16),
                 textposition = grid_polar@axis_labs$adjust,
                 hoverinfo = 'none',
                 showlegend = FALSE,
                 inherit = FALSE) %>%
        # label radial axis
        add_text(x = grid_polar@text_coords$x,
                 y =  -grid_polar@text_coords$y,
                 text = grid_polar@text_coords$text,
                 textposition = 'top center',
                 textfont = list(size = 10),
                 color = I("black"),
                 hoverinfo = 'none',
                 showlegend = FALSE,
                 inherit = FALSE)  %>%

        layout(showlegend = TRUE,
               xaxis = list(title = "", range = c(-grid_polar@r, grid_polar@r),
                            zeroline = FALSE, showline = FALSE,
                            showticklabels = FALSE, showgrid = FALSE,
                            scaleratio = 1, scaleanchor = "y",
                            autoscale = TRUE),
               yaxis = list(title = "", range = c(-grid_polar@r, grid_polar@r),
                            zeroline = FALSE,
                            showline = FALSE,	showticklabels = FALSE,
                            showgrid = FALSE),
               plot_bgcolor = "rgba(0,0,0,0)",
               paper_bgcolor = 'rgba(0,0,0,0)',
               autosize = TRUE)

    return(list("polar"=polar, "cylindrical"=cyl))
}
