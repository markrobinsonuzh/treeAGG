#' replace assays in treeSummarizedExperiment
#'
#' \code{updateTSE} allows to replaces the matrix-like element in assays of
#' \code{treeSummarizedExperiment} with matrix-like element that have the same
#' number of rows but different number of columns.
#'
#' The new \code{assays} have exactly the same rows as the original one that was
#' replaced. But the columns are allowed to be different. The \code{colData} in
#' the original \code{treeSummarizedExperiment} object is removed from the new
#' one.
#'
#' @param assays A list of matrix-like elements.  It should have the same number
#'   of matrix-like elements as the \code{assays} of data specified in the
#'   argument \code{tse}. Also, its matrix-like elements should have  exactly
#'   the same number of rows and rows representing the same thing as those in
#'   the old \code{assays} of \code{tse}. However, the columns could be
#'   different.
#' @param tse A treeSummarizedExperiment object.
#'
#' @export
#' @return A treeSummarizedExperiment object.
#'
#' @author Ruizhu HUANG
#'
updateTSE <- function(assays, tse){

    # rownames
    mat <- assays[[1]]
    rn <- rownames(mat)


    # link data
    linkD <- linkData(tse)

    # keep rows that have rownames included in rn
    # If there are duplicated value in the column nodeLab, use nodeLab_allias
    # to match; otherwise use nodeLab to match.
    if (any(duplicated(linkD$nodeLab))) {
        out_tse <- tse[linkD$nodeLab_allias %in% rn, ]
        ind <- linkData(out_tse)$nodeLab_allias
    } else {
        out_tse <- tse[linkD$nodeLab %in% rn, ]
        ind <- linkData(out_tse)$nodeLab

    }


    # store the output of DA test in the assays.
    mat <- mat[ind, ]
    rownames(mat) <- NULL
    newAssay <- list(mat)

    # update tse
    new_tse <- treeSummarizedExperiment(assays = newAssay,
                                        tree = treeData(out_tse),
                                        rowData = rowData(out_tse),
                                        linkData = linkData(out_tse))
    return(new_tse)
}
