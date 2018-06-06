context("fdr")

test_that("fdr: compute correctly", {
  data("tinyTree")
  # leaf & node level
  expect_equal(fdr(tree = tinyTree, truth = 15,
                   found = 12, level = "leaf"), c(fdr =3/9))
  expect_equal(fdr(tree = tinyTree, truth = 15,
                   found = 12, level = "node"), c(fdr = 6/17))
  expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                   found = c(15, 10), level = "leaf",
                   direction = FALSE), c(fdr = 2/7))
  # direction matters
  expect_equal(fdr(tree = tinyTree, truth = list(16, 13),
                   found = list(13, 16), level = "leaf",
                   direction = TRUE), c(fdr = 1))
  expect_equal(fdr(tree = tinyTree, truth = list(16, 13),
                   found = list(16, 13), level = "leaf",
                   direction = TRUE), c(fdr = 0))

  # node number & node label
  expect_equal(fdr(tree = tinyTree, truth = "Node_15",
                   found = "Node_12", level = "leaf"),
               c(fdr =3/9))
  expect_error(fdr(tree = tinyTree, truth = c("Node_16", "Node_13"),
                   found = c("Node_15", "Node_10"), level = "leaf",
                   direction = FALSE))
  expect_equal(fdr(tree = tinyTree, truth = c("Node_16", "Node_13"),
                   found = c("Node_15", "t3"), level = "leaf",
                   direction = FALSE), c(fdr = 2/7) )
})
