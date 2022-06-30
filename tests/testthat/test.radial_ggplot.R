context("radial_ggplot examples")
library(volcano3D)

test_that("radial_ggplot tests", {
  data(example_data)
  syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
                            data = t(syn_example_rld))
  
  expect_equal(class(radial_ggplot(polar=syn_polar, label_rows=c("COBL"))), 
               c("gg", "ggplot"))
    
})






