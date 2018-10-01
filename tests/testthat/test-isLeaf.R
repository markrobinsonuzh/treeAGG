context("isLeaf")

test_that("shareNode could find correct information", {
    data("tinyTree")
    expect_equal(isLeaf(tree = tinyTree, input = c(5, 4, 18)),
                 c(TRUE, TRUE, FALSE))
    expect_equal(isLeaf(tree = tinyTree, input = c(5, 4, 18)),
                 isLeaf(tree = tinyTree, input = c("t4", "t9", "Node_18" )))
    expect_error(isLeaf(tree = tinyTree, input = 20))
})
