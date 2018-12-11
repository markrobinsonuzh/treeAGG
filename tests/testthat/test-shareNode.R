context("shareNode")

test_that("shareNode could find correct information", {
    data("tinyTree")
    expect_equal(shareNode(node = c('t4','t9'),
                           tree = tinyTree,
                           use.alias = FALSE),
                 c("Node_18" = 18))
    expect_equal(shareNode(node = c('t4','t9', 't1'),
                           tree = tinyTree,
                           use.alias = FALSE),
                 c("Node_16" = 16))
    expect_error(shareNode(node = c('t4','t9', 't11'),
                           tree = tinyTree,
                           use.alias = FALSE))
})
