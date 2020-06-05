context("boxplot_trio examples")
library(volcano3D)
data(example_data)

test_that("boxplot_trio tests", {

   syn_polar <- polar_coords(sampledata = syn_example_meta,
                           contrast = "Pathotype",
                           groups = NULL,
                           pvalues = syn_example_p,
                           expression = syn_example_rld,
                           p_col_suffix = "pvalue",
                           padj_col_suffix = "padj",
                           non_sig_name = "Not Significant",
                           multi_group_prefix = "LRT",
                           significance_cutoff = 0.01,
                           fc_col_suffix='log2FoldChange',
                           fc_cutoff = 0.3)

    t <- boxplot_trio(syn_polar, value = "SLAMF6", plot_method="ggplot")
    expect_equal(length(t), 9)
    
})




