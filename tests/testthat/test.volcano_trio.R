context("volcano_trio examples")
library(volcano3D)
library(volcano3Ddata)

test_that("volcano_treio tests", {
    syn_p_obj <- create_dep(sampledata = syn_metadata, 
               contrast = "Pathotype", 
               pvalues = syn_pvalues,
               p_col_suffix ="pvalue", 
               fc_col_suffix = "log2FoldChange",
               multi_group_prefix = "LRT", 
               expression = syn_rld)
    
    test_obj <- volcano_trio(dep = syn_p_obj,)
    
    expect_equal(length(test_obj), 4)
    expect_equal(names(test_obj), c(rep("", 3), "All"))
    expect_equal(all(unlist(lapply(test_obj, function(x) {
        class(x)[1:2] == c("gg", "ggplot")
        }))), TRUE)
    
})




