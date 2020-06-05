context("volcano3D function examples")
library(volcano3D)
library(volcano3Ddata)

test_that("volcano3D example works", {
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
    
    text_obj <- volcano3D(syn_polar, 
     label_rows = c("FMOD", "LAMP5", "TNNT3"), 
     xy_aspectratio = 1, 
     label_size = 10, 
     z_aspectratio = 0.9)


    expect_equal(class(test_obj)[1], "plotly")
    
    expect_error(volcano3D(iris))
})



