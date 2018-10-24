context("tpr")

test_that("tpr: compute correctly", {
  data("tinyTree")
  # leaf & node level
  expect_equal(tpr(tree = tinyTree, truth = 12,
                   found = c(13, 16), only.Tip = FALSE),
               c(tpr =14/17))
  expect_equal(tpr(tree = tinyTree, truth = 12,
                   found = c(13, 16), only.Tip = TRUE),
               c(tpr = 8/9))
  expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                   found = c(17, 14), only.Tip = TRUE,
                   direction = FALSE), c(tpr = 5/8))

  # direction matters
  expect_equal(tpr(tree = tinyTree, truth = list(16, 13),
                   found = list(17, 14), only.Tip = TRUE,
                   direction = TRUE), c(tpr = 5/8))
  expect_equal(tpr(tree = tinyTree, truth = list(16, 13),
                   found = list(14, 17), only.Tip = TRUE,
                   direction = TRUE), c(tpr = 0))
  expect_equal(tpr(tree = tinyTree, truth = list(16, 13),
                   found = list(c(17, 14), c(NULL)),
                   only.Tip = TRUE,
                   direction = TRUE), c(tpr = 3/8))
  expect_equal(tpr(tree = tinyTree, truth = list(17, 13),
                   found = list(18, 13), only.Tip = FALSE,
                   direction = TRUE), c(tpr = 8/10))

  # node label & node number
  expect_equal(tpr(tree = tinyTree, truth = "Node_12",
                   found = c("Node_13", "Node_16"),
                   only.Tip = FALSE),
               c(tpr =14/17))

})
