#' Table of the top differential nodes
#'
#' \code{topNodes} extracts the top differential nodes in a data frame or a list
#' of data frame for a given pair of group, ranked by p-value or absolute change
#' (e.g. log-fold change)
#'
#' @param data A tree summarizedExperiment object output from
#'   \code{\link{treeAGG}}
#' @param sort.by The name of the column to sort by
#' @param decreasing A logical value. TRUE or FALSE. Should the sorting be
#'   decreasing or not.
#' @param col.rowData The names of columns to be extracted from \code{rowData}.
#' @param col.linkData The names of columns to be extracted from \code{linkData}.
#' @param use.assays A numeric vector. It specifies which table to show. If
#'   NULL, the table correspoinding to the first number storing in the
#'   \code{use.assays} (in \code{metadata} of \code{input} data) will be shown.
#'   To recall, numbers in the \code{use.assays} record which matrix-like elements
#'   in the \code{assaysTable} have been used to do data analysis.
#'
#' @export
#' @return a list of data frame
#' @author Ruizhu HUANG
#'
topNodes <- function(data, sort.by = "FDR", decreasing = FALSE,
                     col.rowData = NULL, col.linkData = NULL,
                     use.assays = NULL){


    if (is.null(use.assays)) {
        use.assays <- 1
    }

    # extract row data
    rData <- rowData(data, internal = TRUE)

    # the name of column that stores results we are interested
    cNam <- paste("result_assay", use.assays, sep = "")
    fData <- rData[, cNam, drop = FALSE]


    if (ncol(fData) == 0) {
        stop("Results can't be found in the internal part of row data")
        }
    # convert to a list: results from different assays as different elements
    list1 <- lapply(seq_len(ncol(fData)), FUN = function(x) {
        fData[, x]

    })
    names(list1) <- colnames(fData)

    # convert to a list: results from different contrasts as different
    # sub-elements
    list2 <- lapply(list1, FUN = function(x) {
        xx <- lapply(seq_len(ncol(x)), FUN = function(y) { x[, y]})
        names(xx) <- colnames(x)
        return(xx)
    })

    # add columns from rowData and linkData
    list3 <- rapply(list2, function(x) {
        xx <- cbind(x, rowData(data)[, col.rowData],
                    linkData(data)[, col.linkData])
        colnames(xx) <- c(colnames(x), col.rowData, col.linkData)
        return(xx)
    }, how = "list")

    # sort according to a column specified by sort.by
    list4 <- rapply(list3, function(x) {
        if (! sort.by %in% colnames(x)) {
            cat("Available columns: ", colnames(x), ".\n" )
            stop("\n ", sort.by,
                 " can't be found. Choose another from the available columns. \n" )
        }
         oo <- order(x[[sort.by]], decreasing = decreasing, na.last = TRUE)
        fx <- x[oo, ]
         return(fx)
    }, how = "list")

    return(list4)


}
