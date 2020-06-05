context("radial_plotly examples")
library(volcano3D)
library(volcano3Ddata)

test_that("radial_plotly tests", {
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
    text_obj <- radial_plotly(polar = syn_polar, label_rows = c("SLAMF6"))
    
    expect_error(radial_plotly(iris))
    expect_equal(unlist(class(test_obj)), c("plotly", "htmlwidget"))
    
})




