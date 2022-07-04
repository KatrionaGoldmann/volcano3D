context("radial_plotly examples")
library(volcano3D)

test_that("radial_plotly tests", {
  data(example_data)
  syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
                            data = t(syn_example_rld))
  
  test_obj <- radial_plotly(polar = syn_polar, label_rows = c("COBL"))
  
  expect_error(radial_plotly(iris))
  expect_equal(unlist(class(test_obj)), c("plotly", "htmlwidget"))
    
})




