context("shareNode")

test_that("shareNode could find correct information", {
    data("tinyTree")
    expect_equal(shareNode(node = c('t4','t9'), tree = tinyTree,
                           return = "label", use.alias = FALSE), "Node_18")
    expect_equal(shareNode(node = c('t4','t9', 't1'), tree = tinyTree,
                           return = "label", use.alias = FALSE), c("Node_16"))
    expect_error(shareNode(node = c('t4','t9', 't11'), tree = tinyTree,
                           return = "label", use.alias = FALSE))
})
