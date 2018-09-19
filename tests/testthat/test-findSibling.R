context("findSibling")

test_that("findSibling could find correct nodes", {
    data(exTree)
    expect_equal(findSibling(tree = exTree, input = 75), 95)
    expect_error(findSibling(tree = exTree, input = 51))
    expect_equal(findSibling(tree = exTree, input = c(75, 76)), c(95, 85))
    expect_equal(findSibling(tree = exTree, input = "t4"), 39)
})
