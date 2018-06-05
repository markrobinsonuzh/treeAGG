context("transNode")

test_that("give errors when the provided arguments are not in correct form", {
  data("tinyTree")

  expect_error(transNode(tree = 2, input = 2))
  expect_error(transNode(tree = tinyTree, input = 20))
})

test_that("transNode could return correct results", {
  data("tinyTree")

  expect_equal(transNode(tree = tinyTree, input = 17), "Node_17")
  expect_equal(transNode(tree = tinyTree, input = "Node_16"), c(Node_16 = 16))
  expect_equal(transNode(tree = tinyTree, input = c(2, 4, 5)), c("t7", "t9", "t4"))

})


