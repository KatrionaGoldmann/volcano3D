#' Three-way radial comparison Polar Plot (using plotly)
#'
#' This function creates an interactive plotly object which maps differential
#' expression onto a polar coordinates.
#' 
#' @param polar A 'volc3d' object with the p-values between groups of interest
#'   and polar coordinates created by \code{\link{polar_coords}},
#'   \code{\link{deseq_polar}} or \code{\link{voom_polar}}.
#' @param type Numeric value whether to use scaled (Z-score) or unscaled (fold
#'   change) as magnitude. Options are 1 = Z-score (default) or 2 =
#'   unscaled/fold change.
#' @param colours A vector of colour names or hex triplets for the
#'   non-significant points and each of the six groups.
#' @param label_rows A vector of row names or numbers to label.
#' @param arrow_length The length of label arrows (default = 50).
#' @param label_size Font size of labels/annotations (default = 14)
#' @param colour_code_labels Logical whether label annotations should be colour
#' coded. If FALSE label_colour is used.
#' @param label_colour HTML colour of annotation labels if not colour coded. 
#' @param grid_colour The colour of the grid (default="grey80")
#' @param grid_width The width of the grid lines (default=1)
#' @param marker_size Size of the markers (default = 6)
#' @param marker_alpha Opacity for the markers (default = 0.7)
#' @param marker_outline_colour Colour for marker outline (default = white)
#' @param marker_outline_width Width for marker outline (default = 0.5)
#' @param axis_title_size Font size for axis titles (default = 16)
#' @param axis_label_size Font size for axis labels (default = 10)
#' @param axis_colour The colour of the grid axes and labels (default="black")
#' @param axis_width The width of the axis lines (default=2)
#' @param axis_ticks A numerical vector of radial axis tick breaks. If
#' NULL this will be calculated using \code{\link[base]{pretty}}.
#' @param axis_angle Angle in radians for the radial axis (default = 5/6).
#' @param ... Optional parameters passed to \code{\link[volcano3D]{polar_grid}}
#' @return Returns a plotly plot featuring variables on a tri-axis radial graph
#' @seealso \code{\link{polar_coords}}
#' @importFrom plotly plot_ly add_trace add_text add_markers layout
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust setNames
#' @importFrom grDevices hsv
#' @importFrom methods is
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://pubmed.ncbi.nlm.nih.gov/31461658/}{
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
#' radial_plotly(polar = syn_polar, label_rows = c("COBL"))


radial_plotly <- function(polar,
                          type = 1,
                          colours = polar@scheme,
                          label_rows = NULL,
                          arrow_length = 50,
                          label_size = 14,
                          colour_code_labels = FALSE,
                          label_colour = "black",
                          grid_colour = "grey80", 
                          grid_width = 1,
                          marker_size = 7,
                          marker_alpha = 0.8,
                          marker_outline_colour = "white",
                          marker_outline_width = 0.5,
                          axis_title_size = 16,
                          axis_label_size = 10,
                          axis_colour = "black",
                          axis_width = 2,
                          axis_ticks = NULL,
                          axis_angle = 5/6,
                          ...){
  if (is(polar, "polar")) {
    args <- as.list(match.call())[-1]
    return(do.call(radial_plotly_v1, args))  # for back compatibility
  }
  if(! is(polar, "volc3d")) stop("Not a 'volc3d' class object")
  df <- polar@df[[type]]
  grid <- polar_grid(df$r, df$z)
  grid@polar_grid <- grid@polar_grid[grid@polar_grid$area != "cylinder", ]
  polar_grid <- grid@polar_grid
  axes <- grid@axes
  axis_labs <- grid@axis_labs
  r <- grid@r
  text_coords <- grid@text_coords
  
  # Annotate gene labels
  if (length(label_rows) != 0) {
    if(! all(is.numeric(label_rows))) {
      if(! all(label_rows %in% rownames(df))) {
        stop("label_rows must be in rownames(df)")
      }}
    if(all(is.numeric(label_rows))) {
      if(! all(label_rows < nrow(df))) {
        stop("label_rows not in 1:nrow(df)")
      }}
    annot <- lapply(label_rows, function(i) {
      row  <- df[i, ]
      theta <- atan(row$y/row$x)
      if(colour_code_labels) ac <- row$col else ac <- label_colour 
      list(x = row$x,
           y = row$y,
           text = rownames(row),
           textangle = 0,
           ax = sign(row$x)*arrow_length*grid@r*cos(theta),
           ay  = -1*sign(row$x)*arrow_length*grid@r*sin(theta),
           font = list(color = ac, size = label_size),
           arrowcolor = ac,
           arrowwidth = 1,
           arrowhead = 0,
           xanchor = ifelse(row$x < 0, "right", "left"),
           xanchor = ifelse(row$y < 0, "bottom", "top"),
           arrowsize = 1.5)
    })
  } else {annot <- list()}
  
  # plotly plot
  # add the grid
  p <- plot_ly(data = df, x = ~x, mode = "none", type = "scatter",
               colors = colours, showlegend = FALSE) %>%
    add_trace(x = polar_grid$x, y = polar_grid$y, color = I(grid_colour),
              line = list(width = grid_width), showlegend = FALSE, 
              type = "scatter", mode = "lines", hoverinfo = "none") %>%
    # add the "horizontal" axes
    add_trace(x = axes$x, y = axes$y, color = I(axis_colour),
              line = list(width = axis_width), showlegend = FALSE, 
              type = "scatter",
              mode = "lines", hoverinfo = "none", inherit = FALSE) %>%
    # add the label text
    add_text(x = axis_labs$x, y = axis_labs$y, 
             text = levels(polar@outcome),
             color = I(axis_colour), type = "scatter", mode = "text",
             textfont = list(size = axis_title_size),
             textposition = axis_labs$adjust, hoverinfo = 'none',
             showlegend = FALSE, inherit = FALSE) %>%
    # label radial axis
    add_text(x = text_coords$x, y = -text_coords$y,
             text = text_coords$text, textposition = 'top center',
             textfont = list(size = axis_label_size), 
             color = I(axis_colour), hoverinfo = 'none', 
             type = "scatter",
             showlegend = FALSE, inherit = FALSE) %>%
    layout(showlegend = TRUE,
           xaxis = list(title = "", range = c(-r*1.02, r*1.4),
                        zeroline = FALSE, showline = FALSE,
                        showticklabels = FALSE, showgrid = FALSE,
                        scaleratio = 1, scaleanchor = "y",
                        autoscale = TRUE),
           yaxis = list(title = "", range = c(-r, r),
                        zeroline = FALSE,
                        showline = FALSE,	showticklabels = FALSE,
                        showgrid = FALSE),
           plot_bgcolor = "rgba(0,0,0,0)",
           paper_bgcolor = 'rgba(0,0,0,0)',
           autosize = TRUE, uirevision=list(editable=TRUE),
           annotations = annot) %>%
    # add the markers
    add_markers(data = df, x = ~x, y = ~y,
                text = rownames(df),
                opacity = marker_alpha,
                color = ~lab,
                colors = colours,
                marker = list(size = marker_size,
                              line = list(color = marker_outline_colour, 
                                          width = marker_outline_width)),
                hoverinfo = 'text', key = rownames(df),
                inherit = FALSE, type = "scatter",
                mode = "markers")
  
  return(p)
  
}
