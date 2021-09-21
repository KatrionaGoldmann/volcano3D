News
=====

# volcano3D 1.3.1
###### 20/09/2021
* Enable cases where one comparison may be a substring of another
* Allow instances with no significant features


# volcano3D 1.3.0
###### 27/03/2021
* Add volcano4D function for rotating/spinning volcano plots
* Add axis and grid width parameters (grid\_width and axis\_width respectively) to radial\_plotly, radial\_ggplot and volcano3D
* Add option (axes\_from\_origin) to allow axes to start from either the origin (default) or the first radial break in polar\_grid. 


# volcano3D 1.2.0
###### 25/02/2021
* Add argument for scene camera in volcano3D

# volcano3D 1.1.0
###### 04/02/2021
* allow colour coding to be based on pvalue or adjusted pvalue according to cutoff_criteria
* allow subsetting by significance groups with the significance_subset function

# volcano3D 1.0.3
###### 15/08/2020
* fix legend dropping levels in volcano_trio
* update volcano3D xy limits to prevent titles disappearing
* allow flexible colour-coding in standard volcano plot

# volcano3D 1.0.2
###### 03/07/2020

* add offset for 3D axis titles
* Allow changes to the hover text 


# volcano3D 1.0.1
###### 26/06/2020

* update vignette for CRAN (remove WebGL)
* update default plotly parameters
* Allow optional colour coding of labels in radial plots
* remove ggplot warnings of NA in geom_path 
* Allow colour of grids and axes to be changed
* Convert 3D labels to annotations 


# volcano3D 1.0.0
###### 29/05/2020

* Combined create\_dep and polar\_coord functions so no longer backwards compatible. 
* Combined the 3D and 2D functions of polar\_grid
* Moved colour selection to individual plotting functions to make it more intuitive
* Made the fold change columns optional
* Allowed custom grid to be passed in
* Pass label columns as plotly keys in volcano3d and radial_plotly
* Pass label columns as plotly keys in volcano3d and radial_plotly
* Add plotly height and width parameters
* Add plotly boxplots
* Create default colours 
* Improve polar_coords speed

# volcano3D 0.1.0
###### 27/04/2020

* This is the initial build of volcano3D
