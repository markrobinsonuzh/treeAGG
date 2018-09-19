context("matTree")

test_that("matTree could find correct path", {
    data(exTree)
    mat <- matTree(tree = exTree)
    expect_equal(dim(mat), c(50, 12))
    expect_equal(mat[1, ], c(1, 55:51, rep(NA, 6)))
})
