context("volcano3D function examples")
library(volcano3D)
library(volcano3Ddata)

test_that("volcano3D example works", {
    syn_p_obj <- create_dep(sampledata = syn_metadata, 
                           contrast = "Pathotype", 
                           pvalues = syn_pvalues,
                           p_col_suffix ="pvalue", 
                           fc_col_suffix = "log2FoldChange",
                           multi_group_prefix = "LRT", 
                           expression = syn_rld)
    
    syn_polar <- polar_coords(dep = syn_p_obj)
    
    test_obj <- volcano3D(syn_polar, 
                         label_size = 10, 
                         xy_aspectratio = 1, 
                         z_aspectratio = 0.9)
    
    
    expect_equal(all(c("SLAMF6", "PARP16") %in% unlist(volcano3D(syn_polar, 
                                   label_rows=c("SLAMF6", "PARP16"),
                                   label_size = 10, 
                                   xy_aspectratio = 1, 
                                   z_aspectratio = 0.9)$x$attrs)), TRUE)

    expect_equal(class(test_obj)[1], "plotly")
    
    expect_error(volcano3D(iris))
})



