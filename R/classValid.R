#' valid leafSummarizedExperiment class
#'
#' \code{validLSE} is to valid treeSummarizedExperiment object.
#'
#'
#' @param object A leafSummarizedExperiment object
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData
#' @importFrom S4Vectors metadata
#' @importFrom utils head
#' @return TRUE or a character string.
#' @keywords internal
#'

# validLSE <- function(object){
#
#     # -----------------------------------------------------------------------
#     # it must have table in assays
#     if (length(assays(object)) < 0) {
#         return("\n there is nothing in assays. \n")
#     }
#
#     # -----------------------------------------------------------------------
#     #  Tree should be a phylo object
#     if (!inherits(metadata(object)$tree, "phylo")) {
#         return("\n tree is not a phylo object")
#     }
#
#     # -----------------------------------------------------------------------
#     # Different leaf nodes are not allowed to use the same labels.
#     tipLab <- metadata(object)$tree$tip.label
#     isDp <- duplicated(tipLab)
#     anyDp <- any(isDp)
#     if (anyDp) {
#         msg <- cat("\n Different leaf nodes using the same label: ",
#                    head(tipLab[isDp])," \n")
#         return(msg)
#     }
#
#     # -----------------------------------------------------------------------
#     # if nodeLab column exist, they should match with the labels of tree
#     leaves
#     nodeLab <- rowData(object)$nodeLab
#     if (!is.null(nodeLab)) {
#         notIn <- any(!nodeLab %in% tipLab)
#         if (notIn) {
#             msg <- cat("\n ", head(setdiff(nodeLab, tipLab)),
#                        " can not be found as labels of tree leaves. \n",
#                    "Check nodeLab column in rowData again.\n")
#             return(msg)
#         }
#     }
#
#     # ------------------------------------------------------------------------
#     # if nodeLab column doesn't exist, rownames should match with the labels
#     # of tree leaves
#     if (is.null(nodeLab)) {
#         rowNam <- rownames(object)
#         notIn <- any(! rowNam %in% tipLab)
#         if (notIn) {
#             msg <- cat("\n ", head(setdiff(rowNam, tipLab)),
#                        " can not be found as labels of tree leaves.
#                    Check rownames again.\n")
#             return(msg)
#         }
#     }
#
#     # -----------------------------------------------------------------------
#     # Note : duplicated value in nodeLab column is allowed because we might
#     # have multiple rows corresponding to a same leaf.
#     return(TRUE)
# }

checkLSE <- function(object){

    errors <- character()
    # -------------------------------------------------------------------------
    # it must have table in assays
    if (length(assays(object)) < 0) {
      msg <- cat("\n there is nothing in assays. \n")
      errors <- c(errors, msg)
    }

    # -------------------------------------------------------------------------
    #  Tree should be a phylo object
    if (!inherits(metadata(object)$tree, "phylo")) {
        msg <- cat("\n tree is not a phylo object")
    }

    # -------------------------------------------------------------------------
    # Different leaf nodes are not allowed to use the same labels.
    tipLab <- metadata(object)$tree$tip.label
    isDp <- duplicated(tipLab)
    anyDp <- any(isDp)
    if (anyDp) {
        msg <- cat("\n Different leaf nodes using the same label: ",
                   head(tipLab[isDp])," \n")
        errors <- c(errors, msg)
    }

    # -------------------------------------------------------------------------
    # if nodeLab column exist, they should match with the labels of tree leaves
    nodeLab <- rowData(object)$nodeLab
    if (!is.null(nodeLab)) {
        notIn <- any(!nodeLab %in% tipLab)
        if (notIn) {
            msg <- cat("\n ", head(setdiff(nodeLab, tipLab)),
                       " can not be found as labels of tree leaves. \n",
                       "Check nodeLab column in rowData again.\n")
            errors <- c(errors, msg)
        }
    }

    # -------------------------------------------------------------------------
    # if nodeLab column doesn't exist, rownames should match with the labels of
    # tree leaves
    if (is.null(nodeLab)) {
        rowNam <- rownames(object)
        notIn <- any(! rowNam %in% tipLab)
        if (notIn) {
            msg <- cat("\n ", head(setdiff(rowNam, tipLab)),
                       " can not be found as labels of tree leaves.
                   Check rownames again.\n")
            errors <- c(errors, msg)
        }
    }

    # -------------------------------------------------------------------------
    # Note : duplicated value in nodeLab column is allowed because we might
    # have multiple rows corresponding to a same leaf.
    if (length(errors) == 0) {TRUE} else {errors}
}
