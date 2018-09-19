context("pathTree")

test_that("pathTree could find correct path", {
    data(tinyTree)
    paths <- pathTree(tinyTree)
    mat <- matTree(tree = tinyTree)
    expect_equal(paths[[1]], c(1, 13:11))

    mm <- mat[1, ]
    expect_equal(paths[[1]], mm[!is.na(mm)])
})
