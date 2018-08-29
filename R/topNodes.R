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
#' @param show.which A numeric vector. It specifies which table to show. If
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
                     show.which = NULL){

    use.assays <- metadata(data)$use.assays
    if (is.null(show.which)) {
        show.which <- use.assays[[1]]
    }
    aData <- assays(data)
    cData <- lapply(seq_along(aData), FUN = function(x) {
        if (x %in% show.which) {

            xx <- cbind(aData[[x]], rowData(data)[ , col.rowData],
                        linkData(data)[, col.linkData])
            colnames(xx) <- c(colnames(aData[[x]]), col.rowData, col.linkData)
            xx
        } else { NULL }

    })

    kData <- lapply(seq_along(cData), FUN = function(x) {
        #y <- x[x$aggKeep, ]
        if (x %in% show.which) {
            y <- cData[[x]]
            oo <- order(y[, sort.by], decreasing = FALSE)
            sy <- y[oo, ]
            return(sy)
        } else { NULL }

    })
    return(kData)


}
