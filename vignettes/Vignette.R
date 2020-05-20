## ----setup, include = FALSE, echo = FALSE-------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
warning = FALSE, 
message = FALSE, 
fig.height = 7, 
fig.width=7, 
fig.align = "center")
library(knitr)
library(kableExtra)

## ---- echo=FALSE--------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(plotly)
library(usethis)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("volcano3D")

## ---- eval = FALSE------------------------------------------------------------
#  library(devtools)
#  install_github("KatrionaGoldmann/volcano3D")

## -----------------------------------------------------------------------------
library(volcano3D)

## ---- echo=FALSE--------------------------------------------------------------
mytable = data.frame(
  sampledata = c("sampledata\ 
  \n\n(required)", 
  "This shows information for each sample in rows and must contain:\ 
                \n * an ID column: Containing the sample IDs. This must be titled ‘ID’.\ 
                \n * a contrast column: A column containing the three-level factor used for contrasts.\ 
                \n \n"), 
  contrast = c("contrast\ 
               \n\n(required)", 
               "The column name in sampledata which contains the three-level factor used for contrast"), 
  pvalues = c("pvalues\ 
              \n\n(required)", 
              "the pvalues data.frame which contains the statistical\
                significance of probes between groups. This contains: \
              \n * three pvalue columns: one for each comparison with \
              column names of format `paste(groups[i], groups[j], p_col_suffix, sep='_')`.\ 
              We recommend \
              using 'limma' or 'DESeq' pipelines to calculate these pvalues for \
              gene expression.\ 
              \n * _optional_ fold change columns: one for each comparison with\
              column names of format `paste0(groups[i], groups[j], fc_col_suffix, sep='_')`  \
              \n * _optional_ adjusted pvalue columns: one for each comparison\
              with column names of format `paste0(groups[i], groups[j], padj_col_suffix, sep='_')` \
              \n * an _optional_ multi-group pvalue column: from a multi-group\
              test with column name of the form `paste0(multi_group_prefix, '_', p_col_suffix)`.This is typically\
              generated using ANOVA or likelihood ratio tests between all three groups. \
              \n * an _optional_ multi-group adjusted pvalue column: from a\
              multi-group test (column names of form `paste0(multiGroupPrefix, '_', padjColSuffix)`). \
              \n For more information on how to create pvalues data frames see the\
              [pvalue generator vignette](https://katrionagoldmann.github.io/volcano3D/articles/pvalues_generator.html)."),
  exp = c("expression\ 
              \n\n(required)", 
          "A data frame or matrix containing the expression data. This \
           is used to calculate z-score and fold change, therefore it should be a \
           normalised expression object such as log transformed or variance stabilised counts."),
  groups = c("groups", "The groups to be compared (in order). If NULL this \
             defaults to `levels(sampledata[, 'contrasts'])`."), 
  pcsuff = c("p_col_suffix", "The suffix of column names with pvalues (default is 'pvalue')."),
  padjsuff = c("padj_col_suffix", "The suffix of column names with adjusted pvalues (default\
  is 'padj'). If NULL the adjusted pvalue is calculated using `p_col_suffix` and\
  `pvalue_method`."),
  padjmeth = c("padjust_method", "The method to calculate adjusted pvalues if not already\
  provided. Must be one of c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH',\
  'BY', 'fdr', 'none'). Default is 'BH'."),
  gc_suff = c("fc_col_suffix", "The suffix of column names with log(fold change) values\
  (default is 'logFC')."),
  mpref = c("multi_group_prefix", "The prefix for columns containing statistics for a\
  multi-group test (this is typically a likelihood ratio test or ANOVA). Default\
  is NULL."),
  lab = c("label_column", "A column name in pvalues which is to be used to label markers\
  of interest at plotting stage. If NULL the rownames will be used.")
)

kable(t(mytable), row.names = FALSE, col.names = c("Variable", "Details")) %>%
  kable_styling(font_size=11)

## -----------------------------------------------------------------------------
data("example_data")

syn_polar <- polar_coords(sampledata = syn_example_meta,
                          contrast = "Pathotype",
                          pvalues = syn_example_p,
                          expression = syn_example_rld,
                          p_col_suffix = "pvalue",
                          padj_col_suffix = "padj",
                          fc_col_suffix = "log2FoldChange",
                          multi_group_prefix = "LRT",
                          non_sig_name = "Not Significant",
                          significance_cutoff = 0.01,
                          label_column = NULL,
                          fc_cutoff = 0.1)

## ---- eval=FALSE--------------------------------------------------------------
#  head(syn_polar@pvalues)

## ---- echo=FALSE--------------------------------------------------------------
head(syn_polar@pvalues) %>%
  kable(col.names = gsub("_l", " _l", gsub("_p", " _p", colnames(syn_polar@pvalues))) )%>%
  kable_styling(font_size=6.5, full_width = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  setNames(data.frame(table(syn_polar@polar$sig)), c("Significance", "Frequency"))

## ---- echo=FALSE--------------------------------------------------------------
table(syn_polar@polar$sig) %>%
  kable(col.names = c("Significance", "Frequency")) %>%
  kable_styling(full_width = F)

## ---- fig.height=2.8, fig.width=7---------------------------------------------
syn_plots <- volcano_trio(polar = syn_polar,
                          sig_names = c("not significant","significant",
                                        "not significant","significant"),
                          colours = rep(c("grey60",  "slateblue1"), 2),
                          text_size = 9,
                          marker_size=1,
                          shared_legend_size = 0.9,
                          label_rows = c("SLAMF6", "BOC", "FMOD"),
                          fc_line = FALSE, 
                          share_axes = FALSE)

syn_plots$All

## ---- fig.height=9, fig.width=9, out.height="150%"----------------------------
radial_plotly(polar = syn_polar, label_rows = c("SLAMF6", "GREM2", "FMOD")) 

## ---- fig.height=4.5, fig.width=7---------------------------------------------
radial_ggplot(polar = syn_polar,
              label_rows = c("SLAMF6", "FMOD", "GREM2"),
              marker_size = 2.3,
              legend_size = 10) +
  theme(legend.position = "right")

