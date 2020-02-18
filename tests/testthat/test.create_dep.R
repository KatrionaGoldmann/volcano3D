context("create_dep function examples")
library(volcano3D)
library(volcano3Ddata)

test_that("create_dep example works", {
    test_obj <- create_dep(sampledata = syn_metadata, 
                           contrast = "Pathotype", 
                           pvalues = syn_pvalues,
                           p_col_suffix ="pvalue", 
                           fc_col_suffix = "log2FoldChange",
                           multi_group_prefix = "LRT", 
                           expression = syn_rld)
    
    expect_equal(test_obj@contrast, "Pathotype")
    expect_equal(class(test_obj)[1], "dep")
    
    expect_error(create_dep(sampledata = syn_metadata, 
                            contrast = "nonsense", 
                            pvalues = syn_pvalues,
                            p_col_suffix ="pvalue", 
                            fc_col_suffix = "log2FoldChange",
                            multi_group_prefix = "LRT", 
                            expression = syn_rld))
    
    expect_error(create_dep(sampledata = iris, 
                            contrast = "nonsense", 
                            pvalues = syn_pvalues,
                            p_col_suffix ="pvalue", 
                            fc_col_suffix = "log2FoldChange",
                            multi_group_prefix = "LRT", 
                            expression = syn_rld))
})



