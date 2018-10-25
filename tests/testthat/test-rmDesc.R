context("rmDesc")

test_that("rmDesc could find correct path", {
    data(tinyTree)
    expect_equal(rmDesc(node = c(5, 4, 18, 3, 2),
                        tree = tinyTree,
                        use.alias = TRUE),
                 c(Node_18 = 18, Leaf_3 = 3, Leaf_2 = 2))
    expect_equal(rmDesc(node = 1:18, tree = tinyTree,
                        use.alias = FALSE),
                 c(Node_11 = 11))
    expect_error(rmDesc(node = 1:20, tree = tinyTree))
})
