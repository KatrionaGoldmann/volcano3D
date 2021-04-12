#' Add mode bar button to rotate the plot
#'
#' @param p The volcano3D plot
#' @param rotate_icon_path The icon path for rotation
#' @param stop_icon_path The icon path for stop button
#' @param rotate_colour The colour for the rotate button (default='lightgreen')
#' @param stop_colour The colour for the stop button (default='tomato')
#' @importFrom htmlwidgets JS onRender
#' @importFrom plotly config
#' @export

volcano4D <- function(p, 
                      rotate_icon_path=NULL, 
                      stop_icon_path=NULL, 
                      rotate_colour='lightgreen', 
                      stop_colour='tomato') {
  
  if(is.null(rotate_icon_path)){
    rotate_icon_path <- "M19.432,7.157c-0.312-1.113-1.624-1.858-3.496-2.17c0.279,0.331,0.532,0.685,0.754,1.06c1.043,0.299,1.748,0.764,1.911,1.344c0.095,0.335,0.014,0.729-0.24,1.169c-0.274,0.476-0.768,1.007-1.455,1.542c0-0.034,0.005-0.067,0.005-0.101c0-3.816-3.094-6.91-6.91-6.91c-3.816,0-6.91,3.094-6.91,6.91c0,1.169,0.293,2.268,0.805,3.232c-1.366-0.277-2.303-0.805-2.495-1.487c-0.094-0.336-0.013-0.729,0.241-1.169c0.138-0.239,0.35-0.496,0.595-0.756c0.011-0.449,0.055-0.89,0.138-1.317c-1.398,1.144-2.112,2.386-1.805,3.476c0.338,1.205,1.845,1.98,3.968,2.24C5.8,15.854,7.774,16.91,10,16.91c3.389,0,6.201-2.44,6.791-5.659C18.735,9.951,19.795,8.448,19.432,7.157 M10,16.047c-1.651,0-3.147-0.664-4.238-1.738c0.147,0.005,0.295,0.008,0.447,0.008c1.502,0,3.195-0.212,4.941-0.658c1.734-0.443,3.297-1.064,4.595-1.776C14.952,14.299,12.682,16.047,10,16.047 M15.998,10.733c-1.27,0.797-2.973,1.554-5.062,2.088c-1.616,0.414-3.251,0.632-4.727,0.632c-0.427,0-0.827-0.025-1.213-0.061C4.338,12.425,3.954,11.258,3.954,10c0-3.339,2.707-6.046,6.046-6.046c3.34,0,6.047,2.708,6.047,6.046C16.047,10.249,16.027,10.492,15.998,10.733"
  }
  
  if(is.null(stop_icon_path)){
    stop_icon_path <- "M512 512H0V0h512v512z"
  }
  
  button_script <- readLines(system.file("/spin_function/spin_button.js",
                                         package = "volcano3D")) %>% 
    gsub("rotate_icon_path", rotate_icon_path, .)  %>% 
    gsub("stop_icon_path", stop_icon_path, .)  %>% 
    gsub("rotate_colour", rotate_colour, .)  %>% 
    gsub("stop_colour", stop_colour, .)  
  
  rotate_button <- list(
    name = "rotate",
    title = 'Rotate',
    icon= list(
      path = rotate_icon_path,
      transform = 'scale(0.84) translate(-1, 0)'
    ),
    val=FALSE, 
    attr='rotating',
    click = htmlwidgets::JS(button_script))
  
  p %>% 
    config(modeBarButtonsToAdd = list(rotate_button)) %>%
    htmlwidgets::onRender(readLines(system.file("spin_function/rotation.js", 
                                                package = "volcano3D")))
}

