context("findOS")

test_that("findOS return error when the ancestor is not in the tree", {
  data("tinyTree")
  expect_error(findOS(ancestor = 20, tree = tinyTree, only.Tip = TRUE))
  expect_error(findOS(ancestor = 2, tree = 2, only.Tip = TRUE))
})

test_that("findOS could find offspring correctly", {
  data("tinyTree")

  expect_setequal(findOS(ancestor = 15, tree = tinyTree,
                         only.Tip = TRUE,
                         self.include = TRUE),
                  c(t9 = 4, t4 = 5, t8 = 6, t10 = 7,
                    t1 = 8, t5 = 9))
  expect_setequal(findOS(ancestor = 15, tree = tinyTree,
                         only.Tip = FALSE, self.include = TRUE),
                  c(t9 = 4, Node_17 = 17, Node_18 = 18,
                    Node_16 = 16, Node_15 = 15, t4 = 5,
                    t8 = 6, t10 = 7, Node_19 = 19, t1 = 8,
                    t5 = 9))

  expect_setequal(findOS(ancestor = 17, tree = tinyTree,
                         only.Tip = FALSE, self.include = FALSE),
                  c(t9 = 4, Node_18 = 18, t4 = 5, t8 = 6))

  expect_setequal(findOS(ancestor = "Node_18", tree = tinyTree,
                         only.Tip = FALSE, self.include = TRUE,
                         use.alias = TRUE),
                  c(Leaf_4 = 4, Node_18 = 18, Leaf_5 = 5))

  expect_setequal(findOS(ancestor = "Node_18", tree = tinyTree,
                         only.Tip = FALSE, self.include = TRUE),
                  c("t9" = 4, "Node_18" = 18, "t4" = 5))

  data("exTree")
  expect_setequal(findOS(ancestor = 94, tree = exTree,
                         only.Tip = TRUE),
                  c(t50 = 43, t10 = 44))

})
