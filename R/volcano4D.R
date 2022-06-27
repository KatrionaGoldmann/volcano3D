#' Add mode bar button to rotate the plot
#'
#' @param p The volcano3D plot
#' @param rotate_icon_path The svg icon path for rotation. If NULL a play 
#' button is used
#' @param stop_icon_path The svg icon path for stop button. If NULL a pause 
#' button is used
#' @param rotate_colour The colour for the rotate button (default="#c7c7c7")
#' @param stop_colour The colour for the stop button (default='#ff6347', 
#' a.k.a 'tomato')
#' @param scale Scaling for rotation button
#' @param speed The rotation speed (default=180) 
#' @param shiny_event_names If using shiny, pass in any shiny event names which
#' should stop rotation when triggered (e.g. shiny_event_names = c('replot'))
#' @importFrom htmlwidgets JS onRender
#' @importFrom plotly config
#' @return Returns a rotating cylindrical 3D plotly plot featuring variables on 
#' a tri-axis radial graph with the -log10(multi-group test p-value) on the 
#' z-axis
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and 
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
#' @importFrom plotly plot_ly add_trace add_text layout %>%
#' @importFrom grDevices hsv rgb
#' @keywords hplot iplot
#' @concept volcanoplot 
#' @export
#' @examples
#' data(example_data)
#' syn_polar <- polar_coords(sampledata = syn_example_meta,
#'                           contrast = "Pathotype", 
#'                           groups = NULL, 
#'                           pvalues = syn_example_p, 
#'                           expression = syn_example_rld, 
#'                           p_col_suffix = "pvalue", 
#'                           padj_col_suffix = "padj", 
#'                           non_sig_name = "Not Significant", 
#'                           multi_group_prefix = "LRT",
#'                           significance_cutoff = 0.01, 
#'                           fc_col_suffix='log2FoldChange',
#'                           fc_cutoff = 0.3)
#' 
#' p <- volcano3D(syn_polar, 
#'     label_rows = c("FMOD", "LAMP5", "TNNT3"), 
#'     xy_aspectratio = 1, 
#'     label_size = 10, 
#'     z_aspectratio = 0.9)
#'     
#' volcano4D(p)
#' @export

volcano4D <- function(p, 
                      rotate_icon_path=NULL, 
                      stop_icon_path=NULL, 
                      rotate_colour="#c7c7c7", 
                      stop_colour='#ff6347', 
                      scale='scale(0.4) translate(-4, -4)', 
                      speed = 180, 
                      shiny_event_names=c()) {
  
  . <- NULL # to appease the CRAN note
  
  if(!inherits(p, "plotly")) stop("Not a plotly plot")
  
  if(is.null(rotate_icon_path)){
    rotate_icon_path <- paste0("M20 33l12-9-12-9v18zm4-29C12.95 4 4 12.95 4", 
                               " 24s8.95 20 20 20 20-8.95 20-20S35.05 4 24 ", 
                               "4zm0 36c-8.82 0-16-7.18-16-16S15.18 8 24", 
                               " 8s16 7.18 16 16-7.18 16-16 16z")
    }
  
  if(is.null(stop_icon_path)){
    stop_icon_path <- paste0("M24 4C12.95 4 4 12.95 4 24s8.95 20 20 20 ", 
                             "20-8.95 20-20S35.05 4 24 4zm-2 28h-4V16h4v16zm8", 
                             " 0h-4V16h4v16z")
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
      transform = scale
    ),
    val=FALSE, 
    attr='rotating',
    click = htmlwidgets::JS(button_script))
  
  js_rotation <- readLines(system.file("spin_function/rotation.js", 
                        package = "volcano3D")) %>%
    gsub('rotation_speed', speed, .) %>%
    gsub("rotate_colour", rotate_colour, .) 
  
  if(length(shiny_event_names) > 0){
    print("catching")
    event_name_str <- paste0("'", paste(shiny_event_names, collapse="', '"), 
                             "'")
    print(event_name_str)
    js_rotation <- gsub("shiny_event_names", event_name_str, js_rotation)
  }
  
  p %>% 
    config(modeBarButtonsToAdd = list(rotate_button)) %>%
    htmlwidgets::onRender(js_rotation)
}

