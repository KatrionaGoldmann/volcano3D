context("volcano_trio examples")
library(volcano3D)
library(volcano3Ddata)

test_that("volcano_trio tests", {
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
                              fc_col_suffix='log2FoldChange',
                              fc_cutoff = 0.3)
    syn_volcano_plots <- volcano_trio(polar=syn_polar)
    expect_equal(class(syn_volcano_plots), "list")
})




