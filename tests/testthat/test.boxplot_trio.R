context("boxplot_trio examples")
library(volcano3D)

test_that("boxplot_treio tests", {
    syn_p_obj <- create_dep(sampledata = syn_metadata, 
               contrast = "Pathotype", 
               pvalues = syn_pvalues,
               p_col_suffix ="pvalue", 
               fc_col_suffix = "log2FoldChange",
               multi_group_prefix = "LRT", 
               expression = syn_rld)
    
    test_obj <- boxplot_trio(dep = syn_p_obj,value="SLAMF6")
    
    expect_equal(length(test_obj), 9)
    expect_equal(class(test_obj)[2], "ggplot")
    
})




