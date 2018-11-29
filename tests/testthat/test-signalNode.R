context("signalNode")

test_that("signalNode could find correct information", {
    data("tinyTree")
    expect_equal(signalNode(node = c('t4','t9'),
                            tree = tinyTree,
                            use.alias = FALSE),
                 c("Node_18" = 18))
    expect_equal(signalNode(node = c('t4','t9', 't1'),
                            tree = tinyTree,
                            use.alias = FALSE),
                 c("Node_18" = 18, "t1" = 8))
    expect_equal(signalNode(node = c('t4','t9', 't1'),
                            tree = tinyTree,
                            use.alias = TRUE),
                 c("Node_18" = 18, "Leaf_8" = 8))
    expect_error(signalNode(node = c('t4','t9', 't11'),
                            tree = tinyTree,
                            use.alias = FALSE))
})
