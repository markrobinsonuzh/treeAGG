context("selNode")

test_that("selNode could find correct information", {
    set.seed(1)
    data("tinyTree")
    toyTable <- matrix(rnbinom(40, size = 1, mu = 10), nrow = 10)
    colnames(toyTable) <- paste(rep(LETTERS[1:2], each = 2), rep(1:2, 2), sep = "_")
    rownames(toyTable) <- tinyTree$tip.label


    dat <- parEstimate(data = toyTable)
    sN <- selNode(tree = tinyTree, data = dat, all = TRUE)
    expect_equal(dim(sN), c(9, 4))
    expect_equal(dim(
        selNode(tree = tinyTree, data = dat, minTip = 4, maxTip = 9,
                minPr = 0, maxPr = 0.8, all = TRUE)), c(2, 4))
    expect_error(selNode(tree = tinyTree, data = dat, minTip = 4, maxTip = 9,
                         minPr = 0, maxPr = 0.2, all = TRUE))
})
