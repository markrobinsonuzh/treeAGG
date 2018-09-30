
# when data is a data frame or a matrix
nodeValue.A <- function(data, fun = sum, tree, message = FALSE) {
    if (!(inherits(data, "data.frame") |
          inherits(data, "matrix"))) {
        stop("data should be a matrix or data.frame")
    }

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (!setequal(rownames(data), tree$tip.label)) {
        chx <- setdiff(rownames(data), tree$tip.label)
        chs <- head(chx)
        stop(cat("The rownames of data don't match the tree tip labels:",
                 chs, "\n"))
    }
    # rename data with the alias of the node label
    rn <- rownames(data)
    leafNum <- transNode(tree = tree, input = rn, message = FALSE)
    leafLab_alias <- transNode(tree = tree, input = leafNum,
                               use.alias = TRUE, message = FALSE)
    rownames(data) <- leafLab_alias

    # nodes
    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)

    ## the node labels
    nN <- length(nodeI)
    nodeLab <- transNode(tree = tree, input = nodeI,
                      use.alias = FALSE,
                      message = FALSE)
    nodeLab_alias <- transNode(tree = tree, input = nodeI,
                         use.alias = TRUE,
                         message = FALSE)

    # calculate counts at nodes
    cNode <- matrix(NA, nrow = nN, ncol = ncol(data))
    rownames(cNode) <- nodeLab_alias

    for (i in seq_len(nN)) {
        node.i <- nodeI[i]
        tips.i <- findOS(ancestor = node.i, tree = tree,
                         only.Tip = TRUE, self.include = TRUE,
                         return = "label", use.alias = TRUE)
        cNode[i, ] <- apply(data[tips.i, ], 2, fun)

        # print out the running process
        if (message) {
           # Sys.sleep(0.001)
            message(i, " out of ", nN , " finished", "\r", appendLF = FALSE)
            # message(i, "\r", appendLF = FALSE)
            flush.console()
        }
    }
    colnames(cNode) <- colnames(data)

    final <- rbind(data, cNode)

    # if there are duplicated value the node label, use the alias of the node
    # labels as the row names
    if (anyDuplicated(c(rn, nodeLab))) {
        rownames(final) <- c(leafLab_alias, nodeLab_alias)
    } else {
        rownames(final) <- c(rn, nodeLab)

    }

    # output
    return(final)
}


# when data is a leafSummarizedExperiment
nodeValue.B <- function(data, fun = sum, message = FALSE) {
    if (!inherits(data, "leafSummarizedExperiment")) {
        stop("\n data should be a leafSummarizedExperiment. \n")
    }


    # extract table and tree from treeSummarizedExperiment.
    tree <- metadata(data)$tree
    table <- assays(data)
    rData <- rowData(data, use.names = FALSE)
    mData <- metadata(data)
    cData <- colData(data)

    # the row names of table
    tipLab <- rData$nodeLab
    if (is.null(tipLab)) { tipLab <- rownames(data)}
    tipNum <- transNode(tree = tree, input = tipLab, message = FALSE)
    tipLab_alias <- transNode(tree = tree, input = tipNum,
                              use.alias = TRUE, message = FALSE)

    # leaves and internal nodes
    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)
    nodeA <- c(leaf, nodeI)

    ## all nodes
    nN <- length(nodeA)
    nodeLab <- transNode(tree = tree, input = nodeA,
                         use.alias = FALSE,
                         message = FALSE)
    nodeLab_alias <- transNode(tree = tree, input = nodeA,
                      use.alias = TRUE,
                      message = FALSE)

    # -------------------------------------------------------------------
    # create empty table nodes
    tableA <- lapply(table, FUN = function(x){
        y <- matrix(0, nrow = nN, ncol = ncol(x))
        colnames(y) <- colnames(data)
        return(y)
    })

    # create rowData for nodes
    rD <- rData[rep(1, nN), ]
    nc <- ncol(rData)

    # use the alias of the node label as the rownames
    rownames(rData) <- tipLab_alias

    # calculate values at nodes
    for (i in seq_len(nN)) {
        node.i <- nodeA[i]
        tips.i <- findOS(ancestor = node.i, tree = tree,
                         only.Tip = TRUE, self.include = TRUE,
                         return = "label", use.alias = TRUE)

        row.i <- match(tips.i, tipLab_alias)

        # if multiple rows of the descendants exist, do calculation as follows
        if (length(row.i) > 1){
            # rowdata: if all rows have the value, keep value; otherwise, use NA
            rdata.i <- rData[row.i, ]
            for (k in seq_len(nc)) {
                iu <- unique(rdata.i[, k])
                rD[i, k] <- ifelse(length(iu) == 1, iu, NA)
                }

            # calculate values (e.g., abundance or intensity) for each node
            for(j in seq_along(tableA)){
                table.j <- table[[j]]
                tableA[[j]][i, ] <- apply(table.j[row.i, , drop = FALSE], 2, fun)
            }
        }

        if (length(row.i) == 1){
            # rowdata
            rD[i, ] <- rData[row.i, ]

            for(j in seq_along(tableA)){
                table.j <- table[[j]]
                tableA[[j]][i, ] <- as.matrix(table.j[row.i, ])
            }
        }

        # print out the running process
        if (message) {
            #Sys.sleep(0.001)
            message(i, " out of ", nN , " finished", "\r", appendLF = FALSE)
            flush.console()
        }
    }

    # update rowdata; column nodeLab is removed
     rdataA <- rD[, !colnames(rD) %in% c("nodeLab")]
     rownames(rdataA) <- NULL

     if (anyDuplicated(nodeLab)) {
         linkD <- DataFrame(nodeLab = nodeLab,
                            nodeLab_alias = nodeLab_alias,
                            nodeNum = nodeA,
                            isLeaf = nodeA %in% leaf,
                            rowID = seq_len(nN))
     } else {
         linkD <- DataFrame(nodeLab = nodeLab,
                            nodeNum = nodeA,
                            isLeaf = nodeA %in% leaf,
                            rowID = seq_len(nN))
     }


    mData.new <- mData[names(mData) != "tree"]
    tse <- treeSummarizedExperiment(linkData = linkD, tree = tree,
                                    assays = tableA, metadata = mData.new,
                                    colData = cData, rowData = rdataA)
    return(tse)

}




#' Calculate entity values at internal nodes
#'
#' \code{nodeValue} calculates value, e.g., count, for each internal
#' node based on the values of its descendant leaves. For example, if the
#' abundance of a cell cluster, which corresponds to an internal node in the
#' tree, is of interest, we could sum the abundance of all cell types in that
#' cluster. Each cell type, in this case, corresponds to a leaf node in the
#' tree. Different calculations, instead of only sum, might be required. This
#' could be achieved by specifying different functions in the argument
#' \strong{fun}.
#'
#' @param data A leafSummarizedExperiment object or a data frame/matrix.
#' @param fun A function to create the value of an internal node based on values
#' at its descendant leaf nodes. The default is sum.
#' @param tree An optional argument. Only need if \code{data} is a data frame or
#'   matrix.
#' @param message A logical value. The default is TRUE. If TRUE, it will print
#'   out the currenet status of a process.
#' @importFrom utils head flush.console
#' @importFrom S4Vectors metadata DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @return A treeSummarizedExperiment if input data is a
#' treeSummarizedExperiment object. A matrix if input data is a matrix or
#' data frame.
#'
#' @export
#' @docType methods
#' @rdname nodeValue
#' @seealso \code{\link{treeSummarizedExperiment}} for treeSummarizedExperiment
#' creation.
#' @author Ruizhu Huang
#' @examples
#' data("tinyTree")
#' if (require(ggtree)) {
#' (p <- ggtree(tinyTree) + geom_text(aes(label = label)))

#' set.seed(1)
#' count <- matrix(rpois(100, 10), nrow =10)
#' rownames(count) <- tinyTree$tip.label
#' colnames(count) <- paste(c("C1_", "C2_"),
#' c(1:5, 1:5), sep = "")
#'
#' count_tinyTree <- nodeValue(data = count,
#' tree = tinyTree, fun = sum)
#'
#' # check to see whether the count of an internal node is the sum
#' # of counts of its descendant leaves.
#' # here, check the first sample as an example
#'
#' nod <- transNode(tree = tinyTree, input = rownames(count_tinyTree),
#'                  message = FALSE)
#' d <- cbind.data.frame(node = nod, count = count_tinyTree[, 1])
#'
#' ggtree(tinyTree) %<+% d + geom_text2(aes(label = count))
#'}
#'
setGeneric("nodeValue", function(data, fun = sum, tree, message = FALSE) {
    standardGeneric("nodeValue")
})

#' @rdname nodeValue
#' @importFrom utils flush.console
setMethod("nodeValue",
          signature(data = "matrix"),
          nodeValue.A)

#' @rdname nodeValue
#' @importFrom utils flush.console
setMethod("nodeValue",
          signature(data = "data.frame"),
          nodeValue.A)

#' @rdname nodeValue
#' @importFrom utils flush.console
#' @importClassesFrom S4Vectors DataFrame
setMethod("nodeValue",
          signature(data = "DataFrame"),
          nodeValue.A)

#' @rdname nodeValue
#' @importFrom utils flush.console
setMethod("nodeValue",
          signature(data = "leafSummarizedExperiment"),
          nodeValue.B)

