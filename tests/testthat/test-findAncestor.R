context("findAncestor")

test_that("findAncestor could find correct nodes", {
    data(exTree)
    expect_equal(findAncestor(tree = exTree, node = c(53, 61),
                              level = 1, use.alias = TRUE),
                 c("Node_52" = 52, "Node_60" = 60))
    expect_equal(findAncestor(tree = exTree, node = 75,
                              level = 2, use.alias = TRUE),
                 c("Node_68" = 68))
    expect_error(findAncestor(tree = exTree, node = 51, level = 1))
    expect_error(findAncestor(tree = exTree, node = 68, level = 3))
    expect_error(findAncestor(tree = exTree, node = 51, level = -1))
})
