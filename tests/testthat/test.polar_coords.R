context("polar_coords examples")
library(volcano3D)


test_that("polar_coords tests", {
    data(example_data)
    syn_polar <- polar_coords(sampledata=syn_example_meta,
                              contrast="Pathotype",
                              groups = NULL,
                              pvalues = syn_example_p,
                              expression = syn_example_rld,
                              p_col_suffix = "pvalue",
                              padj_col_suffix = "padj",
                              fc_col_suffix = NULL,
                              padjust_method = "BH",
                              multi_group_prefix = NULL,
                              non_sig_name = "Not Significant",
                              significance_cutoff = 0.01, 
                              fc_cutoff=0.3, 
                              label_column = NULL)
    
    expect_error(polar_coords(iris))
    expect_equal(class(syn_polar)[1], "polar")
    
})




