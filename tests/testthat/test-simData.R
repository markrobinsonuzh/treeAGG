context("simData")

test_that("doTable could work correctly", {
  data("throat.otu.tab")
  data("throat.tree")

  expect_error(simData(tree = throat.tree,
                       data = throat.otu.tab,
                       ratio = 2))

})
