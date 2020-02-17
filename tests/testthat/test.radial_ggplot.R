context("radial_ggplot examples")
library(volcano3D)

test_that("radial_ggplot tests", {
    syn_p_obj <- create_dep(sampledata = syn_metadata, 
                            contrast = "Pathotype", 
                            pvalues = syn_pvalues,
                            p_col_suffix ="pvalue", 
                            fc_col_suffix = "log2FoldChange",
                            multi_group_prefix = "LRT", 
                            expression = syn_rld)
    
    syn_polar <- polar_coords(dep = syn_p_obj)
    test_obj <- radial_ggplot(polar=syn_polar, fc_cutoff=0.1)
    
    expect_error(radial_plotly(iris))
    expect_equal(unlist(class(test_obj)), c("gg", "ggplot"))
    
})






