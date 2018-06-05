context("findOS")

test_that("findOS return error when the ancestor is not in the tree", {
  data("tinyTree")
  expect_error(findOS(ancestor = 20, tree = tinyTree, only.Tip = TRUE))
  expect_error(findOS(ancestor = 2, tree = 2, only.Tip = TRUE))
})

test_that("findOS could find offspring correctly", {
  data("tinyTree")

  expect_setequal(findOS(
    ancestor = 15,
    tree = tinyTree,
    only.Tip = TRUE,
    self.include = TRUE
  ), 4:9)
  expect_setequal(findOS(
    ancestor = 15,
    tree = tinyTree,
    only.Tip = FALSE,
    self.include = TRUE
  ), c(4, 17, 18, 16, 15, 5, 6, 7, 19, 8, 9))

  expect_setequal(
    findOS(
      ancestor = 17,
      tree = tinyTree,
      only.Tip = FALSE,
      self.include = FALSE
    ),
    c(4, 18, 5, 6)
  )

  expect_setequal(
    findOS(
      ancestor = "Node_18",
      tree = tinyTree,
      only.Tip = FALSE,
      self.include = TRUE,
      return = "node"
    ),
    c(4, 18, 5)
  )
  expect_setequal(
    findOS(
      ancestor = "Node_18",
      tree = tinyTree,
      only.Tip = FALSE,
      self.include = TRUE,
      return = "label"
    ),
    c("t9", "Node_18", "t4")
  )

  data("exTree")
  expect_setequal(findOS(
    ancestor = 94,
    tree = exTree,
    only.Tip = TRUE
  ), c(43:44))

})
