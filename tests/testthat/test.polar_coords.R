context("polar_coords examples")
library(volcano3D)

test_that("polar_coords tests", {
  data(example_data)
  syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
                            data = t(syn_example_rld))
  
  expect_equal(length(unique(syn_polar@df[[1]]$lab) ), 6)
    
})

