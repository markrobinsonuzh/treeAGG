#' Save the analysis results in the rowData of treeSummarizedExperiment object
#'
#' \code{updateTSE} adds the analysis results in the rowData of the
#' treeSummarizedExperiment object. Which tables of the \code{assays} has been
#' used to do data analysis will be specified in the \code{use.assays} of
#' \code{metadata}. The design matrix and contrast vectors are also included in
#' the \code{metadata}.
#'
#' @param result A list with the structure as \code{list(A = list(a1, a2, ...),
#'   B = list(b1, b2, ...))}. The sub-elements, e.g., a1, is a data frame output
#'   from the data analysis. The sub-elements,a1, a2, ..., are obtained by
#'   performing analysis on the same table A in the \code{assays} of \code{tse}
#'   with different contrasts. b1 and a1 are obtained by perfoming analysis on
#'   different tables (A and B, respectively) in the \code{assays} of \code{tse}
#'   with the same contrast.
#' @param tse A treeSummarizedExperiment object.
#' @param use.assays A numeric vector. It specifies which tables are used for
#'   the analysis. The order should match with the order of \code{results} (A
#'   and B in the example structure)
#' @param design A design matrix that is used in the analysis to get
#'   \code{result}.
#' @param contrast A list of contrast vectors that are used in the analysis to
#'   get \code{result}. The order should match with the subelements of the
#'   \code{results}, like a1, a2, ..., in the example structure.
#' @param fit NULL. If using \code{\link{edgeR}} package for the analysis, it
#'   could store the output from \code{\link[edgeR]{glmFit}}. A object of
#'   \code{\link[edgeR]{DGEGLM-class}}
#' @export
#' @return A treeSummarizedExperiment object.
#'
#' @author Ruizhu HUANG
#'
#' @examples
#'
#' library(edgeR)
#' library(S4Vectors)
#' set.seed(1)
#' y <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(y) <- paste(rep(LETTERS[1:3], each = 10), rep(1:10,3), sep = "_")
#' rownames(y) <- tinyTree$tip.label
#'
#' rowInf <- DataFrame(nodeLab = rownames(y),
#'                     var1 = sample(letters[1:3], 10, replace = TRUE),
#'                     var2 = sample(c(TRUE, FALSE), 10, replace = TRUE))
#' colInf <- DataFrame(gg = factor(sample(1:3, 30, replace = TRUE)),
#'                     group = rep(LETTERS[1:3], each = 10))
#' toy_lse <- leafSummarizedExperiment(tree = tinyTree, rowData = rowInf,
#'                                     colData = colInf,
#'                                     assays = list(y, (2*y), 3*y))
#'
#' toy_tse <- nodeValue(data = toy_lse, fun = sum, tree = tinyTree, message = TRUE)
#'
#' # extract the abundace table
#' count <- assays(toy_tse, use.nodeLab = TRUE)[[1]]
#'
#' # The sample size is the sum of cell counts of clusters on the leaf level of the tree.
#' tipCount <- assays(toy_tse[linkData(toy_tse)$isLeaf,], use.nodeLab = TRUE)[[1]]
#' libSize <- apply(tipCount, 2, sum)
#'
#' # create DGEList
#' y <- DGEList(counts = count, lib.size = libSize,
#'              remove.zeros = FALSE)
#'
#' # calculate normalisation factors
#' y <- calcNormFactors(object = y, method = "TMM")
#'
#' # construct design matrix
#' sample_inf <- colData(toy_tse)
#' design <- model.matrix(~ gg + group, data = sample_inf)
#'
#' # estimate dispersion
#' y <- estimateGLMRobustDisp(y, design = design)
#'
#' # fit the negative binomial GLMs
#' fit <- glmFit(y, design = design, prior.count = 0.125)
#'
#' # run likelihood ratio tests
#' # contrast
#' contrast1 <- c(0, 0, 0, -1, 1)
#' contrast2 <- c(0, -1, 1, 0, 0)
#'
#' # contrast is not specified here, so the last coefficient is tested.
#' lrt1 <- glmLRT(fit, contrast = contrast1)
#' lrt2 <- glmLRT(fit, contrast = c (0, 0, 1, 0, 0))
#' # matC <- cbind(c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
#' # lrt3 <- glmLRT(fit, contrast = matC)
#'
#' # extract table
#' tab1 <- topTags(lrt1, n = Inf, adjust.method = "BH", sort.by = "none")$table
#' tab2 <- topTags(lrt2, n = Inf, adjust.method = "BH", sort.by = "none")$table
#' #tab3 <- topTags(lrt3, n = Inf, adjust.method = "BH", sort.by = "none")$table
#'
#' # put the analysis result in the toy_tse
#' contrastList <- list(contrastG1 = contrast1, contrastG2 = contrast2)
#' #,contrastG3 = matC)
#'
#' # the analaysis is performed on the first table of assays
#' # add more elements if there are results from more tables of assays
#' res <- list(assay1 = list(contrastG1 = tab1, contrastG2 = tab2))
#' #, contrastG3 = tab3))
#'
#' new_tse <- updateTSE(result = res, tse = toy_tse,
#' use.assays = 1, design = design, contrast = contrastList,
#' fit = fit)
#'
#' # the results is added to a column called result_assay1
#' rowData(toy_tse)
#' rowData(new_tse)
#'
updateTSE <- function(result, tse, use.assays, design, contrast, fit = NULL){

    # ------------------------------- check -----------------------------------
    # the input result should be a list
    # a list
    #  - A                (result gained from A count table)
    #  -- a1, a2, a3, a4  (each gained from A under contrast c1, c2, c3, c4)
    #  - B                (result gained from B count table)
    #  -- b1, b2, b3, b4  (each gained from B under contrast c1, c2, c3, c4)
    # example: list(A = list(a1, a2, a3, a4), B = list(b1, b2, b3, b4))
    if (!inherits(result, "list")) {
            stop("\n The input result should be a list. \n")
    }
    sub_result <- unlist(result, recursive = FALSE)
    class_subelements <- unlist(lapply(sub_result, class))

    if (any(!class_subelements %in% c("data.frame", "DataFrame"))) {
        stop("\n Each sub-element of the input result should be a data frame. \n")
    }

    # the number of elements should equal to the length of use.assays
    # The example structure has A & B (length 2)
    if (length(result) != length(use.assays)) {
        stop("\n The length of *result* should match the length of *use.assays*. \n")
    }

    # the number of sub-elements should be the same
    # The example structure has a1-a4 (length 4) and b1-b4 (length 4)
    le <- unlist(lapply(result, length))
    if (length(unique(le)) != 1) {
        stop("\n The sub-elements in *result* have different length: ", le, "\n")
    }

    # if the contrast is NULL
    if (is.null(contrast)) {
        contrast <- lapply(seq_along(result[[1]]), function(x) {NULL})
        names(contrast) <- "contrastNULL"
    }

    # the number of sub-elements should equal to the length of contrast
    # The example structure has a1-a4 (length 4)
    if (unique(le) != length(contrast)) {
        stop("\n The number of sub-elements in the *result* (", unique(le),
             ") isn't equal to the length of contrast (",
             length(contrast), "). \n")
    }

    # names are required for the contrast
    if (is.null(names(contrast))) {
        stop("The input contrast should be named. \n")
    }

    # the names of sub-elements should match the contrast
    namList <- lapply(result, names)
    namList <- unlist(namList[!duplicated(namList)])
    check.name <- all(namList == names(contrast))

    if (!check.name) {
        stop("\n The names for the sub-elements of *result* should match with the
             names of the *contrast*. \n")
    }
    # names the result list
    tt1 <- result
    names(tt1) <- paste("result_assay",
                        use.assays,
                        sep = "" )
    # tt1 <- lapply(tt1, function(x) {
    #     names(x) <- names(contrast)
    #     return(x)
    #     })

    # add NA rows for the deleted rows in the analysis
    # This is to make sure the result table has the same rows as the count table
    # in the assays of tse
    tt2 <- lapply(seq_along(tt1), function(x) {
            # add and rearrange rows so that the output is also row-wise
            # corresponding to the assays data.
            res <- lapply(tt1[[x]], function(y) {
                # find rows deleted
                idx <- rownames(y)
                idc <- setdiff(linkData(tse)$nodeLab_alias, idx)

                df <- DataFrame(y)
                dm <- matrix(NA, nrow = length(idc), ncol = ncol(df),
                             dimnames = list(idc, colnames(df)))
                dfc <- DataFrame(dm)
                dfA <- rbind(df, dfc)
                od <- linkData(tse)$nodeLab_alias
                if (is.null(od)) {
                    od <- linkData(tse)$nodeLab
                }
                dfA <- dfA[od,]
                return(dfA)
            })
            return(res)
    })

    # reshape: convert a list into a dataFrame
    tt3 <- lapply(seq_along(tt2), FUN = function(j) {
            x <- tt2[[j]]
            dx <- x[[1]][, 0]
            for (i in seq_along(x)) {
                xi <- x[[i]]
                nxi <- names(x)[i]
                if (is.null(nxi)) { nxi <- "contrastNULL"}
                dx[[nxi]] <- xi

            }
            return(dx)
    })
    names(tt3) <- names(tt1)

    # put the analysis result in the rowData
    # set the class as internal_rowData
    for (i in seq_along(tt3)) {
        if (!is.null(tt3[[i]])) {
            rowData(tse)[[names(tt3)[i]]] <- as(tt3[[i]], "internal_rowData")
        }
    }

    # put the contrast, matrix, use.assays in the metadata
    metadata(tse)$output_glmFit <- fit
    metadata(tse)$use.assays <- use.assays
    metadata(tse)$design <- design
    metadata(tse)$contrast <- contrast

    # output the new tse
    return(tse)

}


