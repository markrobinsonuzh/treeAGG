context("rmDesc")

test_that("rmDesc could find correct path", {
    data(tinyTree)
    aV <- c(5, 4, 18, 3, 2)
    rn <- rmDesc(node = aV, tree = tinyTree)
    expect_equal(rn, c(18, 3, 2))
    expect_equal(rmDesc(node = 1:18, tree = tinyTree), 11)
    expect_error(rmDesc(node = 1:20, tree = tinyTree))
})
