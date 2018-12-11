context("matTree")

test_that("matTree could find correct path", {
    data(exTree)
    mat <- matTree(tree = exTree)
    vv <- c(1, 55:51, rep(NA, 6))
    names(vv) <- paste("L", seq_along(vv), sep = "")
    expect_equal(dim(mat), c(50, 12))
    expect_equal(mat[1, ], vv)
})
