context("polar_coords examples")
library(volcano3D)
library(volcano3Ddata)

test_that("polar_coords tests", {
    syn_p_obj <- create_dep(sampledata = syn_metadata, 
                            contrast = "Pathotype", 
                            pvalues = syn_pvalues,
                            p_col_suffix ="pvalue", 
                            fc_col_suffix = "log2FoldChange",
                            multi_group_prefix = "LRT", 
                            expression = syn_rld)
    
    syn_polar <- polar_coords(dep = syn_p_obj)
    
    expect_error(polar_coords(iris))
    expect_equal(class(syn_polar)[1], "polar")
    
})




