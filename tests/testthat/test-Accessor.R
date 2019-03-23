context("classAccessor")

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
# leafSummarizedExperiment
lse <- leafSummarizedExperiment(assays = list(toyTable),
                                rowData = rowInf,
                                colData = colInf,
                                tree = tinyTree)
# treeSummarizedExperiment
tse <- treeSummarizedExperiment(assays = list(toyTable),
                                rowData = rowInf,
                                colData = colInf,
                                tree = tinyTree)


test_that("assays could extract table successfully", {
              expect_equal(assays(lse)[[1]], assays(tse)[[1]])
              expect_equal(assays(lse)[[1]], toyTable)
              expect_equal(assays(tse)[[1]], toyTable)
              # expect_equal(rownames(assays(tse, use.nodeLab = TRUE,
              #                     withDimnames = TRUE)[[1]]), rowInf$nodeLab)
              expect_equal(rownames(assays(tse, use.nodeLab = FALSE,
                                           withDimnames = TRUE)[[1]]),
                           rownames(toyTable))
          })

test_that("assays could be written successfully", {
              assays(lse)[[2]] <- 2*toyTable
              assays(tse)[[2]] <- 2*toyTable
              expect_equal(assays(lse)[[2]], assays(tse)[[2]])
              expect_equal(assays(lse)[[2]], 2*toyTable)
              expect_equal(assays(tse)[[2]], 2*toyTable)
          })

test_that("row data could be extracted successfully", {
    expect_equal(rowData(lse), rowInf)
    expect_equal(colnames(rowData(tse)), setdiff(colnames(rowInf), "nodeLab"))
    expect_equal(rowData(tse), rowInf[, !colnames(rowInf) %in% "nodeLab"])
    expect_equal(rowData(tse, use.names = TRUE), rowInf[, !colnames(rowInf) %in% "nodeLab"])
    expect_equal(rownames(rowData(tse, use.names = FALSE)), NULL)
})

test_that("row data could be written successfully", {
    lse1 <- lse
    rowData(lse1)$test <- rep(1, nrow(lse))

    tse1 <- tse
    rowData(tse1)$test <- rep(1, nrow(tse))
    expect_equal(rowData(lse1)$test, rep(1, nrow(lse)))
    expect_equal(rowData(tse1)$test, rep(1, nrow(tse)))
})


test_that("column data could be extracted successfully", {
    expect_equal(colData(lse), colInf)
    expect_equal(colData(tse), colInf)
})
