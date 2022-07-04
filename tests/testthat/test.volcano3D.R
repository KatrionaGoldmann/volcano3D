context("volcano3D function examples")
library(volcano3D)

test_that("volcano3D example works", {
  data(example_data)
  syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
                            data = t(syn_example_rld))
  
  test_obj <- volcano3D(syn_polar)
  
  expect_true(inherits(test_obj, "plotly"))
  
  expect_error(volcano3D(iris))
})



