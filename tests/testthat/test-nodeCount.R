context("nodeCount")


test_that("nodeCount requires correct input argument", {

  data("tinyTree")
  count <- matrix(1:100, nrow =10)
  rownames(count) <- tinyTree$tip.label
  colnames(count) <- paste(c("C1_", "C2_"),
  c(1:5, 1:5), sep = "")

  expect_error(nodeCount(tipTable = count, tree = NULL,
                         stree = NULL, fun = sum))
  expect_equal(nodeCount(tipTable = count, tree = tinyTree,
                         stree = NULL, fun = sum)[18, 1], 9)
  expect_equal(nodeCount(tipTable = count, tree = tinyTree,
                         stree = NULL, fun = sum)[nrow(count)+1, 1], sum(count[,1]))
  expect_equal(nodeCount(tipTable = count, tree = tinyTree,
                         stree = NULL, fun = mean)[nrow(count)+1, 1], mean(count[,1]))

})
