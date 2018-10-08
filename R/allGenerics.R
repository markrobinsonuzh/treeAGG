# -----------------------------------------------------------------------------
### Accessors for leafSummarizedExperiment
# -----------------------------------------------------------------------------

#' @importMethodsFrom SummarizedExperiment assayNames
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("assayNames", signature("leafSummarizedExperiment"),
          function(x, ...){
              callNextMethod()
          })

#' @importMethodsFrom SummarizedExperiment "assayNames<-"
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("assayNames<-", signature("leafSummarizedExperiment"),
          function(x, ..., value){
              callNextMethod()
          })

#' @importMethodsFrom SummarizedExperiment assays
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("assays", signature("leafSummarizedExperiment"),
          function(x, ..., withDimnames = TRUE){
              callNextMethod()
          })

#' @importMethodsFrom SummarizedExperiment "assays<-"
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("assays<-", signature("leafSummarizedExperiment"),
          function(x, ..., withDimnames = TRUE, value){
              callNextMethod()
          })

#' @importMethodsFrom SummarizedExperiment rowData
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("rowData", "leafSummarizedExperiment",
          function(x, use.names = TRUE,...) {
    callNextMethod()
})


#' @importMethodsFrom SummarizedExperiment "rowData<-"
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("rowData<-", "leafSummarizedExperiment",
          function(x, ..., value) {
    callNextMethod()
})

#' @importMethodsFrom SummarizedExperiment colData
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("colData", signature("leafSummarizedExperiment"),
          function(x, ...){
              callNextMethod()
          })

#' @importMethodsFrom SummarizedExperiment "colData<-"
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("colData<-", signature(x = "leafSummarizedExperiment"),
          function(x, ..., value){
              callNextMethod()
          })



#' @importMethodsFrom S4Vectors metadata
#' @importFrom methods callNextMethod
#' @rdname leafSummarizedExperiment-accessor
#' @export
setMethod("metadata", signature("leafSummarizedExperiment"),
          function(x) {
              callNextMethod()
          } )


# -----------------------------------------------------------------------------
### Accessors for treeSummarizedExperiment
# -----------------------------------------------------------------------------
#' @importMethodsFrom SummarizedExperiment assayNames
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("assayNames", signature("treeSummarizedExperiment"),
          function(x, ...){
              callNextMethod()
          })

#' @importMethodsFrom SummarizedExperiment "assayNames<-"
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("assayNames<-", signature("treeSummarizedExperiment", "character"),
          function(x, ..., value){
              callNextMethod()
          })


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

#' @importMethodsFrom SummarizedExperiment "assays<-"
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("assays<-", signature(x = "treeSummarizedExperiment"),
          function(x, ..., withDimnames = TRUE, value){
              callNextMethod()
          })

#' @importMethodsFrom SummarizedExperiment rowData
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("rowData", "treeSummarizedExperiment",
          function(x, use.names = TRUE, internal = TRUE, ...) {
    if (internal) {
        callNextMethod()
    } else {
        vv <- callNextMethod()
        cv <- unlist(lapply(vv, class))
        isInternal <- cv == "internal_rowData"
        vv[, !isInternal, drop = FALSE]
    }

})


#' @importMethodsFrom SummarizedExperiment "rowData<-"
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("rowData<-", "treeSummarizedExperiment",
          function(x, ..., value) {
    callNextMethod()
})



#' @importMethodsFrom SummarizedExperiment colData
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("colData", signature("treeSummarizedExperiment"),
          function(x, ...){
              callNextMethod()
          })

#' @importMethodsFrom SummarizedExperiment "colData<-"
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("colData<-", signature(x = "treeSummarizedExperiment"),
          function(x, ..., value){
              callNextMethod()
          })


#' @importMethodsFrom S4Vectors metadata
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("metadata", signature("treeSummarizedExperiment"),
          function(x) {
              callNextMethod()
          } )

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
setMethod("[", signature(x = "treeSummarizedExperiment"),
          function(x, i, j){
                xx <- SummarizedExperiment(assays = assays(x),
                                         rowData = rowData(x),
                                         colData = colData(x),
                                         metadata = metadata(x))
              out <- callNextMethod(x = xx, i = i, j = j)

              linkD <- x@linkData
              if (missing(i)) {
                  i <- seq_len(nrow(linkD))
              }

              linkD <- linkD[i, , drop = FALSE]

              # final <- new("treeSummarizedExperiment", out, linkData = linkD,
              #              tree = x@tree)
              final <- treeSummarizedExperiment(tree = x@treeData,
                                                linkData = linkD,
                                                assays = assays(out),
                                                rowData = rowData(out),
                                                colData = colData(out),
                                                metadata = metadata(out))
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



