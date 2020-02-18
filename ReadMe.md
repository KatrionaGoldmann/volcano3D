[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/volcano3D)](https://cran.r-project.org/package=volcano3D)
[![Codecov test
coverage](https://codecov.io/gh/r-lib/volcano3D/branch/master/graph/badge.svg)](https://codecov.io/gh/r-lib/volcano3D?branch=master)
[![Downloads](https://cranlogs.r-pkg.org/badges/volcano3D?color=blue)](https://cran.rstudio.com/package=volcano3D)

volcano3D <img src="logo.png" align="right" alt="" width="120" />
=================================================================

The volcano3D package enables exploration of probes differentially
expressed between three groups. Its main purpose is for the
visualisation of differentially expressed genes in a three-dimensional
volcano plot. These plots can be converted to interactive visualisations
using plotly.

The vignette explores a case study from the PEAC rheumatoid arthritis
trial (Pathobiology of Early Arthritis Cohort). The methodology has been
published in [‘Lewis, Myles J., et al. “Molecular portraits of early
rheumatoid arthritis identify clinical and treatment response
phenotypes.” Cell reports 28.9 (2019):
2455-2470.’](https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1)
with an interactive web tool available at <https://peac.hpc.qmul.ac.uk>.

This tool acts as a searchable interface to examine relationships
between individual synovial and blood gene transcript levels and
histological, clinical, and radiographic parameters, and clinical
response at 6 months. An interactive interface allows the gene module
analysis to be explored for relationships between modules and clinical
parameters. The PEAC interactive web tool was creating as an [R Shiny
app](https://shiny.rstudio.com) and deployed to the web using a server.

Installing and Loading the package
----------------------------------

#### Prerequisite

-   ggplot2
-   ggpubr
-   ggrepel
-   plotly

### Install from CRAN

*Not yet publicly available:*

    install.packages("volcano3D")

### Install from Github

*Not yet publicly available:*

    library(devtools)
    install_github("KatrionaGoldmann/volcano3D")

------------------------------------------------------------------------

Citation
========

volcano3D was developed by the bioinforamtics team at the [Experimental
Medicine & Rheumatology department](https://www.qmul.ac.uk/whri/emr/)
and [Centre for Translational
Bioinformatics](https://www.qmul.ac.uk/c4tb/) at Queen Mary University
London.

If you use this package please cite as:

> Lewis, Myles J., et al. “Molecular portraits of early rheumatoid
> arthritis identify clinical and treatment response phenotypes.” Cell
> reports 28.9 (2019): 2455-2470.

or using:

    citation("volcano3D")

    ## 
    ## To cite package 'volcano3D' in publications use:
    ## 
    ##   Katriona Goldmann and Myles Lewis (2020). volcano3D: Interactive
    ##   Plotting of Three-Way Differential Expression Analysis. R package
    ##   version 0.1.0.9000. https://github.com/KatrionaGoldmann/volcano3D
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {volcano3D: Interactive Plotting of Three-Way Differential Expression
    ## Analysis},
    ##     author = {Katriona Goldmann and Myles Lewis},
    ##     year = {2020},
    ##     note = {R package version 0.1.0.9000},
    ##     url = {https://github.com/KatrionaGoldmann/volcano3D},
    ##   }
    ## 
    ## ATTENTION: This citation information has been auto-generated from the
    ## package DESCRIPTION file and may need manual editing, see
    ## 'help("citation")'.
