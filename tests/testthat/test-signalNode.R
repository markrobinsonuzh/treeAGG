context("signalNode")

test_that("signalNode could find correct information", {
    data("tinyTree")
    expect_equal(signalNode(node = c('t4','t9'), tree = tinyTree,
                           label = TRUE), "Node_18")
    expect_equal(signalNode(node = c('t4','t9', 't1'), tree = tinyTree,
                           label = TRUE), c("Node_18", "t1"))

    expect_error(signalNode(node = c('t4','t9', 't11'), tree = tinyTree,
                            label = TRUE))
})
