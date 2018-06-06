context("tpr")

test_that("tpr: compute correctly", {
  data("tinyTree")

  expect_equal(tpr(tree = tinyTree, truth = 12,
                   found = c(13, 16), level = "node"), c(tpr =14/17))
  expect_equal(tpr(tree = tinyTree, truth = 12,
                   found = c(13, 16), level = "leaf"), c(tpr = 8/9))
  expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                   found = c(17, 14), level = "leaf",
                   direction = FALSE), c(tpr = 5/8))
  # direction matters
  expect_equal(tpr(tree = tinyTree, truth = list(16, 13),
                   found = list(17, 14), level = "leaf",
                   direction = TRUE), c(tpr = 5/8))
  expect_equal(tpr(tree = tinyTree, truth = list(16, 13),
                   found = list(14, 17), level = "leaf",
                   direction = TRUE), c(tpr = 0))
  expect_equal(tpr(tree = tinyTree, truth = list(16, 13),
                   found = list(c(17, 14), c(NULL)), level = "leaf",
                   direction = TRUE), c(tpr = 3/8))
  expect_equal(tpr(tree = tinyTree, truth = list(17, 13),
                   found = list(18, 13), level = "node",
                   direction = TRUE), c(tpr = 8/10))

})
