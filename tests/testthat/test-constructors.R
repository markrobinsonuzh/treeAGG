context("classConstructor")

# assays data
set.seed(1)
toyTable <- matrix(rnbinom(40, size = 1, mu = 10), nrow = 10)
colnames(toyTable) <- paste(rep(LETTERS[1:2], each = 2), rep(1:2, 2), sep = "_")
rownames(toyTable) <- paste("entity", seq_len(10), sep = "")

data("tinyTree")
# row data
rowInf <- DataFrame(var1 = sample(letters[1:2], 10, replace = TRUE),
                    var2 = sample(c(TRUE, FALSE), 10, replace = TRUE),
                    nodeLab = tinyTree$tip.label,
                    row.names = rownames(toyTable))
# column data
colInf <- DataFrame(gg = c(1, 2, 3, 3),
                    group = rep(LETTERS[1:2], each = 2),
                    row.names = colnames(toyTable))


# treeSummar

test_that("leafSummarizedExperiment constuctor works", {
    # leafSummarizedExperiment
    lse <- leafSummarizedExperiment(assays = list(toyTable),
                                    rowData = rowInf,
                                    colData = colInf,
                                    tree = tinyTree)

    expect_is(lse, "leafSummarizedExperiment")
    expect_is(leafSummarizedExperiment(rowData = rowInf,
                                       colData = colInf,
                                       tree = tinyTree),
              "leafSummarizedExperiment")
})

test_that("treeSummarizedExperiment constuctor works", {
    # leafSummarizedExperiment
    tse <- treeSummarizedExperiment(assays = list(toyTable),
                                    rowData = rowInf,
                                    colData = colInf,
                                    tree = tinyTree)

    expect_is(tse, "treeSummarizedExperiment")
    expect_is(treeSummarizedExperiment(rowData = rowInf,
                                       colData = colInf,
                                       tree = tinyTree),
              "treeSummarizedExperiment")
})
