context("boxplot_trio examples")
library(volcano3D)
data(example_data)

test_that("boxplot_trio tests", {

  syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
                            data = t(syn_example_rld))
  
  trio <- boxplot_trio(syn_polar, value = "COBL", plot_method="ggplot")
  expect_equal(length(trio), 9)
    
})

