context("findAncestor")

test_that("findAncestor could find correct nodes", {
    data(exTree)
    expect_equal(findAncestor(tree = exTree, node = c(53, 61), level = 1),
                 c(52, 60))
    expect_equal(findAncestor(tree = exTree, node = 75, level = 2),
                 68)
    expect_error(findAncestor(tree = exTree, node = 51, level = 1))
    expect_error(findAncestor(tree = exTree, node = 68, level = 3))
    expect_error(findAncestor(tree = exTree, node = 51, level = -1))
})
