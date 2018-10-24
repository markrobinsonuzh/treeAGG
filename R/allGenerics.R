
# -----------------------------------------------------------------------------
### Accessors for treeSummarizedExperiment
# -----------------------------------------------------------------------------
#' @importMethodsFrom SummarizedExperiment assays
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("assays", signature("treeSummarizedExperiment"),
          function(x, use.nodeLab = FALSE, ..., withDimnames = FALSE){
              out <- callNextMethod(x, withDimnames)

              if (use.nodeLab) {
                  nodeLab <- x@linkData$nodeLab
                  if (any(duplicated(nodeLab))) {
                      nodeLab <- x@linkData$nodeLab_alias
                  }

                  outR <- lapply(out, function(x) {
                      rownames(x) <- nodeLab
                      x
                  })
              } else {
                  outR <- out
              }
              return(outR)
          })


#' @importMethodsFrom SummarizedExperiment rowData
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("rowData", "treeSummarizedExperiment",
          function(x, use.names = TRUE, ...) {
              vv <- callNextMethod()
              cv <- unlist(lapply(vv, class))
              isInternal <- cv == "internal_rowData"
              vv[, !isInternal, drop = FALSE]
          })

#' @rdname treeSummarizedExperiment-accessor
#' @export
setGeneric("linkData", function(x) {
    standardGeneric("linkData")
})


#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("linkData", signature("treeSummarizedExperiment"),
          function(x) {
              x@linkData
          })

#' @rdname treeSummarizedExperiment-accessor
#' @export
setGeneric("treeData", function(x) {
    standardGeneric("treeData")
})

#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("treeData", signature("treeSummarizedExperiment"),
          function(x) {
              x@treeData
          })

#' @importFrom methods callNextMethod
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom S4Vectors metadata
#' @rdname treeSummarizedExperiment-accessor
#' @export
# @importFrom BiocGenerics replaceSlots
setMethod("[", signature(x = "treeSummarizedExperiment"),
          function(x, i, j){
              # subset slots inherited from SummarizedExperiment
              xx <- SummarizedExperiment(assays = assays(x),
                                         rowData = rowData(x),
                                         colData = colData(x),
                                         metadata = metadata(x))
              out <- callNextMethod(x = xx, i = i, j = j)

              # subset new slots
              linkD <- x@linkData
              if (missing(i)) {
                  i <- seq_len(nrow(linkD))
              }
              linkD <- linkD[i, , drop = FALSE]

              # final <- BiocGenerics:::replaceSlots(x,
              #                                      linkData = linkD,
              #                                      assays = assays(out),
              #                                      rowData = rowData(out),
              #                                      colData = colData(out),
              #                                      metadata = metadata(out))

              final <- new("treeSummarizedExperiment", out,
                           treeData = x@treeData,
                           linkData = linkD)
              return(final)
          })


#' @keywords internal
#' @importFrom methods callNextMethod
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "treeSummarizedExperiment", function(object) {
    callNextMethod()
    cat(
        "treeData:", " a phylo \n",
        "linkData:", " a ", class(linkData(object)), " with ",
        ncol(linkData(object)), " columns \n",
        sep=""
    )
})




