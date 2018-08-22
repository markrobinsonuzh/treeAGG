
# -----------------------------------------------------------------------------
### Documentation of the accessor function
### Function code is in the file allGenerics.R
# -----------------------------------------------------------------------------




# -----------------------------------------------------------------------------
### leafSummarizedExperiment
# -----------------------------------------------------------------------------

#' Accessor functions for leafSummarizedExperiment
#'
#' Accessor functions to extract different elements from
#' \strong{leafSummarizedExperiment} object.
#'
#' @param x A leafSummarizedExperiment object
#' @param ... For assay, ... may contain withDimnames, which is forwarded to
#'   assays. For rowData, arguments passed through ... are forwarded to mcols.
#'   For cbind, rbind, ... contains leafSummarizedExperiment objects to be
#'   combined. For other accessors, ignored.
#' @param withDimnames A logical(1), indicating whether dimnames should be
#'   applied to extracted assay elements. Setting withDimnames=FALSE increases
#'   the speed and memory efficiency with which assays are extracted.
#'   withDimnames=TRUE in the getter assays<- allows efficient complex
#'   assignments (e.g., updating names of assays, names(assays(x,
#'   withDimnames=FALSE)) = ... is more efficient than names(assays(x)) = ...);
#'   it does not influence actual assignment of dimnames to assays.
#' @param value	An object of a class specified in the S4 method signature. See
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} for more details.
#' @name leafSummarizedExperiment-accessor
#' @return Elements from \code{leafSummarizedExperiment}.
#' @author Ruizhu HUANG
#' @seealso \code{\link{treeSummarizedExperiment}}
#'   \code{\link{treeSummarizedExperiment-accessor}}
#'   \code{\link{leafSummarizedExperiment}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
NULL

# -----------------------------------------------------------------------------
### treeSummarizedExperiment
# -----------------------------------------------------------------------------
#' Accessor functions for treeSummarizedExperiment
#'
#' Accessor functions to extract different elements from
#' \strong{treeSummarizedExperiment} object.
#'
#' @param x A treeSummarizedExperiment object
#' @param ... For assay, ... contains \code{use.nodeLab}, which is forwarded to
#'   assays. For rowData, arguments passed through ... are forwarded to mcols.
#'   For cbind, rbind, ... contains leafSummarizedExperiment objects to be
#'   combined. For other accessors, ignored.
#' @param use.nodeLab A logical(1), indicating whether the rownames of assay
#'   elements should use node labels (the column \code{nodeLab} in
#'   \code{linkData} if there is not duplicated values; otherwise the column
#'   \code{nodeLab_allias} in \code{linkData} is used.)
#' @param withDimnames A logical(1), indicating whether dimnames should be
#'   applied to extracted assay elements. Setting withDimnames=FALSE increases
#'   the speed and memory efficiency with which assays are extracted.
#'   withDimnames=TRUE in the getter assays<- allows efficient complex
#'   assignments (e.g., updating names of assays, names(assays(x,
#'   withDimnames=FALSE)) = ... is more efficient than names(assays(x)) = ...);
#'   it does not influence actual assignment of dimnames to assays.
#' @param value	An object of a class specified in the S4 method signature. See
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} for more details.
#' @param i,j The subscripts that can act to subset the rows and columns of
#'   \code{x}, that is the matrix elements of assays.
#' @name treeSummarizedExperiment-accessor
#' @return Elements from \code{treeSummarizedExperiment}.
#' @author Ruizhu HUANG
#' @seealso \code{\link{treeSummarizedExperiment}}
#'   \code{\link{treeSummarizedExperiment-accessor}}
#'   \code{\link{leafSummarizedExperiment}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
NULL