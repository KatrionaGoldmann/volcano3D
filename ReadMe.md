[![Lifecycle:
Maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![License: GPL
v2](https://img.shields.io/badge/License-GPL%20v2-mediumpurple.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![CRAN
status](https://www.r-pkg.org/badges/version/volcano3D)](https://cran.r-project.org/package=volcano3D)
[![Downloads](https://cranlogs.r-pkg.org/badges/volcano3D?color=orange)](https://cran.rstudio.com/package=volcano3D)
[![2020-06-16](https://img.shields.io/badge/last%20git%20commit-2020--06--16-turquoise.svg)](https://github.com/KatrionaGoldmann/volcano3D/blob/master/NEWS.md)
[![HitCount](http://hits.dwyl.com/KatrionaGoldmann/volcano3D.svg)](http://hits.dwyl.com/KatrionaGoldmann/volcano3D)

volcano3D <img src="logo.png" align="right" alt="" width="200" hspace="20" />
=============================================================================

The volcano3D package enables exploration of probes differentially
expressed between three groups. Its main purpose is for the
visualisation of differentially expressed genes in a three-dimensional
volcano plot. These plots can be converted to interactive visualisations
using plotly.

The
[vignette](file:///Users/kgoldmann/Documents/Analyses/volcano_package/volcano3D/docs/articles/Extended_Vignette.html)
explores a case study from the PEAC rheumatoid arthritis trial
(Pathobiology of Early Arthritis Cohort). The methodology has been
published in [Lewis, Myles J., et al. *Molecular portraits of early
rheumatoid arthritis identify clinical and treatment response
phenotypes*. Cell reports 28.9 (2019): 2455-2470. (DOI:
10.1016/j.celrep.2019.07.091)](https://doi.org/10.1016/j.celrep.2019.07.091)
with an interactive web tool available at <https://peac.hpc.qmul.ac.uk>.

This tool acts as a searchable interface to examine relationships
between individual synovial and blood gene transcript levels and
histological, clinical, and radiographic parameters, and clinical
response at 6 months. An interactive interface allows the gene module
analysis to be explored for relationships between modules and clinical
parameters. The PEAC interactive web tool was creating as an [R Shiny
app](https://shiny.rstudio.com) and deployed to the web using a server.

There are also supplementary vignettes for further information on:

-   [setting up the input pvalues data
    frame](https://katrionagoldmann.github.io/volcano3D/articles/pvalues_generator.html)
-   [using the volcano3D package to create and deploy a shiny
    app](https://katrionagoldmann.github.io/volcano3D/articles/shiny_builder.html)

Getting Started
---------------

### Prerequisites

-   [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
-   [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
-   [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)
-   [plotly](https://cran.r-project.org/web/packages/plotly/index.html)

### Install from CRAN

    install.packages("volcano3D")

### Install from Github

    library(devtools)
    install_github("KatrionaGoldmann/volcano3D")
    library(volcano3D)

### volcano3D data

The sample data can then also be installed either from
[source](https://github.com/KatrionaGoldmann/volcano3Ddata) or using:

    install_github("KatrionaGoldmann/volcano3Ddata")

Citation
--------

volcano3D was developed by the bioinformatics team at the [Experimental
Medicine & Rheumatology department](https://www.qmul.ac.uk/whri/emr/)
and [Centre for Translational
Bioinformatics](https://www.qmul.ac.uk/c4tb/) at Queen Mary University
London.

If you use this package please cite as:

    citation("volcano3D")

    ## 
    ## To cite package 'volcano3D' in publications use:
    ## 
    ##   Katriona Goldmann and Myles Lewis (2020). volcano3D: Interactive
    ##   Plotting of Three-Way Differential Expression Analysis.
    ##   https://katrionagoldmann.github.io/volcano3D/index.html,
    ##   https://github.com/KatrionaGoldmann/volcano3D.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {volcano3D: Interactive Plotting of Three-Way Differential Expression
    ## Analysis},
    ##     author = {Katriona Goldmann and Myles Lewis},
    ##     year = {2020},
    ##     note = {https://katrionagoldmann.github.io/volcano3D/index.html,
    ## https://github.com/KatrionaGoldmann/volcano3D},
    ##   }

or:

> Lewis, Myles J., et al. *Molecular portraits of early rheumatoid
> arthritis identify clinical and treatment response phenotypes*. Cell
> reports 28.9 (2019): 2455-2470.
