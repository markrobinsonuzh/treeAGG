
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
                  # The pair between rows of assays and rows of linkData
                  rowID <- x@linkData$rowID
                  pl <- lapply(seq_along(rowID), FUN = function(x) {
                     xx <- rowID[[x]]
                     cbind(iLink = rep(x, length(xx)), iAssay = xx)
                  })
                  pm <- do.call(rbind, pl)
                  ps <- pm[order(pm[, "iAssay"], decreasing = FALSE), ]
                  ind <- ps[, "iLink"]
                  nodeLab <- x@linkData$nodeLab
                  if (any(duplicated(nodeLab))) {
                      lab <- x@linkData$nodeLab_alias[ind]
                  } else {
                      lab <- nodeLab[ind]
                  }

                  outR <- lapply(out, function(x) {
                      rownames(x) <- lab
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
#  @importFrom BiocGenerics normalize
#' @importFrom S4Vectors metadata
#' @rdname treeSummarizedExperiment-accessor
#' @export
#'
setMethod("[", signature(x = "treeSummarizedExperiment"),
          function(x, i, j){

             # Subset the traditional slots from SummarizedExperiment
              nx <- callNextMethod()

              # Extract the map information between rows of assays and
              # nodes of the tree
              lk <- x@linkData
              rid <- lk$rowID
              nodeNum <- lk$nodeNum
              node <- lapply(seq_along(rid), FUN = function(x) {
                  rep(x, length(rid[[x]]))
              })
              pair <- cbind(rid = unlist(rid), nodeID = nodeNum[unlist(node)])
              pair <- pair[order(pair[, "rid"], decreasing = FALSE), ]

              if (!missing(i)) {
                  pair <- pair[i, ]
              }

             # Update the column 'rowID' in new slot (linkData)
              rowID <- lapply(seq_along(nodeNum), FUN = function(x) {
                  nID <- pair[, "nodeID"]
                  which(nID %in% x)

              })
              lk$rowID <- rowID

              # update slots
              final <- BiocGenerics:::replaceSlots(nx,
                                                   linkData = lk,
                                                   treeData = x@treeData)
              return(final)
          })


#' @rdname treeSummarizedExperiment-accessor
#' @export
setGeneric("getByLink", function(x, subset) {
    standardGeneric("getByLink")
})

#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("getByLink", signature(x = "treeSummarizedExperiment"),
          function(x, subset) {
              # Access the link data
              lx <- x@linkData

              # Get the column 'rowID'
              e <- substitute(subset)
              lID <- base::subset(lx, subset = eval(e),
                                  select = "rowID", drop = TRUE)
              rID <- unlist(lID)
              # # Subset the object
              x[rID, ]
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




