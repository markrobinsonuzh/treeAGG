context("signalNode")

test_that("signalNode could find correct information", {
    data("tinyTree")
    expect_equal(signalNode(node = c('t4','t9'), tree = tinyTree,
                           return = "label", use.alias = FALSE),
                 "Node_18")
    expect_equal(signalNode(node = c('t4','t9', 't1'), tree = tinyTree,
                            return = "label", use.alias = FALSE),
                 c("Node_18", "t1"))
    expect_equal(signalNode(node = c('t4','t9', 't1'), tree = tinyTree,
                            return = "number", use.alias = FALSE),
                 c(18, 8))
    expect_equal(signalNode(node = c('t4','t9', 't1'), tree = tinyTree,
                            return = "number", use.alias = TRUE),
                 c(18, 8))
    expect_error(signalNode(node = c('t4','t9', 't11'), tree = tinyTree,
                            return = "label", use.alias = FALSE))
})
