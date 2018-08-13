#-------------------------------------------------------------------------------
#' A S3 class for the phylogenetic tree
#'
#' The \pkg{ape} package does not export its phylo class, probably because it
#' is not really defined formally anywhere. Technically, it is an S3 class
#' extended from the class list. Any exported definitions from the \pkg{ape}
#' package would be preferred to use if available.
#' @importFrom methods setOldClass
#' @rdname phylo-class
#' @keywords internal
phylo <- structure(list(), class = "phylo")
setOldClass("phylo")

#-------------------------------------------------------------------------------
# treeSummarizedExperiment
#-------------------------------------------------------------------------------
#' An S4 class treeSummarizedExperiment
#'
#' The class \code{treeSummarizedExperiment} is extended from the class
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} by including a
#' \strong{phylo} object in the metadata. It further restricts the rowData to
#' have a column named \strong{nodeLab}. This column provides the labels of
#' nodes on the tree that are corresponding to the rows of matrix in
#' \strong{assays}. If multiple rows of a matrix in \strong{assays} are matched
#' to a same node on the tree, then the same node label is shared by these rows.
#' @importFrom methods setClass
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @name treeSummarizedExperiment-class
#' @rdname treeSummarizedExperiment-class
#' @exportClass treeSummarizedExperiment
#' @include classValid.R
#'
setClass("treeSummarizedExperiment",
         contains = "SummarizedExperiment",
         validity = validTSE)


#-------------------------------------------------------------------------------
#' A virtual class to combine "matrix", "data.frame".
#' @importFrom methods setClassUnion
#' @rdname matrixDataframe
#' @keywords internal
setClassUnion("matrixDataframe", c("matrix", "data.frame"))



