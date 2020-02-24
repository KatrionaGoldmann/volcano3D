setClassUnion("character_or_NULL", c("character", "NULL"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))

#' An S4 class to define the polar grid and coordinates for polar
#' differential expression plots
#' 
#' @slot sampledata Sample data with column ID and contrast
#' @slot contrast The column name in `sampledata`` which contains the
#'   three-group contrast factor used for comparisons.
#' @slot multi_group_test Column name prefix for statistical tests between
#'   all three groups
#' @slot pvalues A data frame containing the p-values, adjusted p-values,
#'   and log2(fold changes) for all three comparisons between groups in the 
#'   contrast factor, as well as optional multi-group tests.
#' @slot expression An optional data frame or matrix containing the
#'   expression data
#' @slot polar A data.frame containing:
#'   \itemize{
#'       \item{The mean expression for each of the three groups in comparison}
#'       \item{The z-score polar coordinates: 'y_zscore', 'x_zscore' and 
#'       'r_zscore'}
#'       \item{The fold-change polar coordinates: 'y_fc', 'x_fc' and 'r_fc'}
#'       \item{'angle': The angle in radians for polar coordinates}
#'       \item{'angle_degrees': The angle in degrees}
#'       \item{'hue': The continuous colour for volcano3D}
#'       \item{'col': The discrete colour for volcano3D}
#'       \item{'maxExp': The maximally expressed group}
#'       \item{'sig': The significance group}
#'   }
#' @slot non_sig_name The category name for variables which are classed as not
#' significant
setClass("polar", slots = list(sampledata = "data.frame",
                               contrast = "character",
                               multi_group_test = "character_or_NULL",
                               pvalues = "data.frame",
                               expression = "df_or_matrix",
                               polar = "df_or_matrix",
                               non_sig_name = "character"))


#' Coordinates for Three Way Polar Plot
#'
#' This function creates a data frame for downstream polar plots containing
#' the p-values from a three-way group comparison. 
#' @param dep A dep object with the pvalues between groups of interest. Created
#' by \code{\link{create_dep}}.
#' @param primary_colours The colours for the primary groups. I.e. points up in
#' only one of the three groups (default = c('blue', 'red', 'green3')).
#' @param secondary_colours The colours for points in secondary groups. I.e.
#' those up in two of the three groups (default = c("violet", "gold3", "cyan")).
#' @param non_sig_colour Colour of non-significant points
#' @param non_sig_name Name to assign non-significant points
#' @param significance_cutoff Value defining the significance cut-off (pvalues
#' below this will be coloured \code{non_sig_colour})
#' @return Returns an S4 polar object containing:
#' \itemize{
#'   \item{'polar'} A data.frame containing:
#'   \itemize{
#'       \item{The mean expression for each of the three groups in comparison}
#'       \item{The z-score polar coordinates: 'y_zscore', 'x_zscore' and 
#'       'r_zscore'}
#'       \item{The fold-change polar coordinates: 'y_fc', 'x_fc' and 'r_fc'}
#'       \item{'angle': The angle in radians for polar coordinates}
#'       \item{'angle_degrees': The angle in degrees}
#'       \item{'hue': The continuous colour for volcano3D}
#'       \item{'col': The discrete colour for volcano3D}
#'       \item{'maxExp': The maximally expressed group}
#'       \item{'sig': The significance group}
#'   }
#'   \item{'pvalues'} A data frame containing the p-values, adjusted p-values,
#'   and log2(fold changes) for all three comparisons between groups in the 
#'   contrast factor, as well as optional multi-group tests.
#'   \item{'sampledata'} Sample data with column ID and contrast
#'   \item{'contrast'} The column name in `sampledata`` which contains the
#'   three-group contrast factor used for comparisons.
#'   \item{'multi_group_test'} Column name prefix for statistical tests between
#'   all three groups
#'   \item{'expression'} An optional data frame or matrix containing the
#'   expression data
#'   \item{'non_sig_name'} The category name for variables which are classed as 
#'   not significant
#' }
#' @keywords pvalue, polar, plot
#' @importFrom grDevices col2rgb hsv
#' @references
#' Lewis, Myles J., et al. (2019).
#' \href{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}{
#' Molecular portraits of early rheumatoid arthritis identify clinical and
#' treatment response phenotypes.}
#' \emph{Cell reports}, \strong{28}:9
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
#' syn_polar <- polar_coords(dep = syn_p_obj)
#' table(syn_polar@polar$sig) 

polar_coords <- function(dep,
                         primary_colours = c('green3', 'blue', 'red'),
                         secondary_colours = c("cyan", "purple", "gold2"),
                         non_sig_colour = "grey60",
                         non_sig_name = "Not Significant",
                         significance_cutoff = 0.01){
    
    expression <- dep@expression
    polar_pvalues <- dep@pvalues
    
    if(! class(expression) %in% c("data.frame", "matrix")) {
        stop("expression must be a data frame or matrix")
    }
    if(! class(polar_pvalues) %in% "data.frame") {
        stop("pvalues must be a data frame")
    }
    if(class(try(col2rgb(non_sig_colour),silent = TRUE)) == "try-error") {
        stop('non_sig_colour must be a valid colour')
        
    }
    if(any(unlist(lapply(primary_colours, function(x) {
        class(try(col2rgb(x), silent = TRUE)) == "try-error"
    })))) stop('all primary_colours must be valid colours')
    if(any(unlist(lapply(secondary_colours, function(x) {
        class(try(col2rgb(x), silent = TRUE)) == "try-error"
    })))) stop('all secondary_colours must be valid colours')
    
    sampledata <- dep@sampledata
    contrast <- dep@contrast
    if(! "ID" %in% colnames(sampledata)) {
        stop("There is no ID column in the sampledata")
    }
    if(length(primary_colours) != 3) {
        stop('The colour vector must be of length 3')
    }
    if(length(secondary_colours) != 3) {
        stop('The colour vector must be of length 3')
    }
    if(! is.numeric(significance_cutoff)) {
        stop('significance_cutoff must be a numeric')
    }
    if(! (significance_cutoff >= 0 & significance_cutoff <= 1)) {
        stop('significance_cutoff must be between 0 and 1')
    }
    
    # Check alignement
    if(identical(as.character(sampledata$ID), colnames(expression)) == FALSE) {
        stop("sampledata and expression data not properly aligned")
    }
    
    # Capture the expression score for each group (fc and z-score)
    expression_scaled <- t(scale(t(expression)))
    polar_colours <- as.data.frame(t(apply(expression_scaled, 1, function(x) {
        tapply(x, droplevels(sampledata[, dep@contrast]), mean, na.rm = TRUE)
    })))
    polar_coloursFC <- as.data.frame(t(apply(expression, 1, function(x) {
        tapply(x, droplevels(sampledata[, dep@contrast]), mean, na.rm = TRUE)
    })))
    
    if(! identical(rownames(expression), rownames(dep@pvalues))) {
        stop('expression and pvalues not aligned (rownames are not identical)')
    }
    if(! identical(rownames(polar_colours), rownames(dep@pvalues))) {
        stop('expression and pvalues not aligned (rownames are not identical)')
    }
    
    polar_colours$Name <- rownames(polar_colours)
    sampledata$contrast <- droplevels(sampledata[, dep@contrast])
    contrast_groups <- levels(sampledata[, dep@contrast])
    if(length(contrast_groups) != 3) {
        stop("The number of variables in the contrast column of sampledata
             does not equal 3")
    }
    
    # Calculate the polar coordinates (uses radians)
    polar_colours$y_zscore <- sinpi(1/3)*(polar_colours[, contrast_groups[2]] -
                                        polar_colours[, contrast_groups[3]])
    polar_colours$x_zscore <- polar_colours[, contrast_groups[1]] -
        (cospi(1/3)*(polar_colours[, contrast_groups[3]] +
                         polar_colours[, contrast_groups[2]]))
    polar_colours$y_fc <- sinpi(1/3)*(polar_coloursFC[,contrast_groups[2]] -
                                          polar_coloursFC[,contrast_groups[3]])
    polar_colours$x_fc <- polar_coloursFC[,contrast_groups[1]] -
        (cospi(1/3)*(polar_coloursFC[,contrast_groups[3]] +
                         polar_coloursFC[,contrast_groups[2]]))
    
    # angle runs from -0.5 to +0.5
    polar_colours$angle <- atan2(polar_colours$y_zscore,
                                 polar_colours$x_zscore)/(2*pi)
    polar_colours$angle_degrees <- (polar_colours$angle %% 1)*360
    
    # rotate and modulus
    polar_colours$angle <- (polar_colours$angle + 2/3) %% 1
    
    # Calculate the colours by significance
    raw_hue <- hsv(polar_colours$angle[!is.na(polar_colours$angle)], 1, 1)
    polar_colours$hue <- 'grey'
    #polar_colours$raw_hue <- within(polar_colours, hue[!is.na(angle)] <- )
    if(! is.null(dep@multi_group_test)){
        polar_colours$hue[polar_pvalues[, paste(dep@multi_group_test, "padj")]
                          >= significance_cutoff] <- non_sig_colour
    }
    polar_colours$r_zscore <- with(polar_colours, sqrt(x_zscore^2 + y_zscore^2))
    polar_colours$r_fc <- with(polar_colours, sqrt(x_fc^2 + y_fc^2))
    
    
    # Set up the colours - pick the most highly expressed group
    polar_colours$col <- as.character(primary_colours[max.col(
        polar_colours[, 1:3])])
    if(! is.null(dep@multi_group_test)){
        polar_colours$col[polar_pvalues[, 
                                        paste(dep@multi_group_test, "padj")] >= 
                              significance_cutoff] <- non_sig_colour
    }
    polar_colours$maxExp <- colnames(polar_colours)[max.col(
        polar_colours[, 1:3])]
    comp_map <- colnames(polar_colours)[1:3]
    comp_cols <- colnames(polar_pvalues)[grepl("pvalue", 
                                               colnames(polar_pvalues))
                                         & colnames(polar_pvalues) != 
                                             paste(dep@multi_group_test, 
                                                   "pvalue")]
    
    polar_colours$col[polar_pvalues[,comp_cols[1]] >= significance_cutoff &
                      polar_pvalues[,comp_cols[2]] >= significance_cutoff &
                      polar_pvalues[,comp_cols[3]] >= significance_cutoff] <-
        non_sig_colour
    polar_colours$maxExp[polar_pvalues[,comp_cols[1]] >= significance_cutoff &
                         polar_pvalues[,comp_cols[2]] >= significance_cutoff &
                         polar_pvalues[,comp_cols[3]] >= significance_cutoff] <-
        non_sig_name
    
    
    
    
    # Calculate which significance group each gene belongs to
    index <- polar_colours$col != non_sig_colour & ! is.na(polar_colours$col)
    if (any(index != FALSE)){
        pairwise_comp <- data.frame(polar_pvalues[index, comp_cols],
                                    row.names = rownames(polar_pvalues)[index])
        colnames(pairwise_comp) <- gsub(" P.Value", "", comp_cols)
        
        # Determine which groups are significant
        pairwise <- data.frame(ifelse(pairwise_comp <= significance_cutoff, 
                                      "1", 
                                      "0"))
        pairwise$min <- as.numeric(
            apply(polar_colours[index, 1:3], 1, function(x) {
                as.numeric(which.min(x))
            }))
        pairwise$min2 <- comp_map[match(colnames(polar_colours)[1:3],
                                        comp_map)][pairwise$min]
        
        # create a string determiing: min (of ABC),  sig AvB,  sig BvC, sig CvA
        pairwiseSig <- paste0(pairwise$min2, pairwise[,1],
                              pairwise[,2], pairwise[,3])
        pairwise$col <- non_sig_colour
        sigRes <- c(paste0(comp_map[2], "1..|", comp_map[3], "..1"),
                    paste0(comp_map[1], "1..|", comp_map[3], ".1."),
                    paste0(comp_map[1], "..1|", comp_map[2], ".1."))
        maxRes <- c(paste0(comp_map[3], ".11"),
                    paste0(comp_map[1], "1.1"),
                    paste0(comp_map[2], "11."))
        pairwise$maxExp <- non_sig_name
        
        # Fill in colours
        for(iter in 1:3){
            pairwise$col[grep(sigRes[iter], pairwiseSig)] <- 
                primary_colours[iter]
            pairwise$maxExp[grep(sigRes[iter], pairwiseSig)] <-
                colnames(polar_colours)[iter]
        }
        
        for(iter in 1:3){
            pairwise$col[grep(maxRes[iter], pairwiseSig)] <- 
                secondary_colours[iter]
            iterIDs <- sort(c(iter, iter%%3 + 1))
            pairwise$maxExp[grep(maxRes[iter], pairwiseSig)] <-
                paste(colnames(polar_colours)[iterIDs[1]], "+",
                      colnames(polar_colours)[iterIDs[2]], "+", sep = "")
        }
        
        pairwise$Name <- rownames(pairwise)
        polar_colours$col <- pairwise$col[match(rownames(polar_colours),
                                                rownames(pairwise))]
        polar_colours$sig <- pairwise$maxExp[match(rownames(polar_colours),
                                                   rownames(pairwise))]
    } else polar_colours$col <- non_sig_colour
    
    polar_colours$col[is.na(polar_colours$col)] <- non_sig_colour
    polar_colours$hue[polar_colours$col == non_sig_colour] <- non_sig_colour
    polar_colours$sig[is.na(polar_colours$sig)] <- non_sig_name
    
    polar_colours$col <- factor(polar_colours$col)
    polar_colours$sig <- factor(polar_colours$sig)
    
    polar_colours <- polar_colours[, c("Name",
                                       comp_map,
                                       "y_zscore", "x_zscore", "r_zscore",
                                       "x_fc", "y_fc", "r_fc",
                                       "angle",
                                       "angle_degrees",
                                       "hue",
                                       "col",
                                       "maxExp",
                                       "sig")]
    
    if(! identical(rownames(polar_colours), rownames(dep@pvalues))) {
        stop('Misalignment of rows in polar and pvalues')
    }
    
    methods::new("polar",
                 polar = polar_colours,
                 pvalues = dep@pvalues,
                 sampledata = dep@sampledata,
                 contrast = dep@contrast,
                 multi_group_test = dep@multi_group_test,
                 expression  = dep@expression,
                 non_sig_name = non_sig_name)
} 






