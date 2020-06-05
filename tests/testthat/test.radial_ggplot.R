context("radial_ggplot examples")
library(volcano3D)
library(volcano3Ddata)

test_that("radial_ggplot tests", {
    data(example_data)
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
                           fc_cutoff = 0.3)
    expect_equal(class(radial_ggplot(polar=syn_polar, label_rows=c("SLAMF6"))), 
                 c("gg", "ggplot"))
    
})






