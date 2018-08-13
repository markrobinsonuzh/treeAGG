

#' Calculate entity values for internal nodes
#'
#' \code{nodeValue} calculates the entity value at each internal node.
#'
#' @param tree A phylo object
#' @param data A matrix or data frame.
#' @param fun A function.
#'
#' @importFrom utils head
#' @keywords internal
#' @return A count table (matrix class) with a row representing a node and a
#'   column representing a sample.
#' @author Ruizhu Huang
#'

nodeValue.A <- function(tree, data, fun = sum) {
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

    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)

    ## calculate counts for nodes
    nN <- length(nodeI)
    nNam <- transNode(tree = tree, input = nodeI)

    # calculate counts at nodes
    cNode <- matrix(NA, nrow = nN, ncol = ncol(data))
    rownames(cNode) <- nNam
    for (i in seq_len(nN)) {
        node.i <- nodeI[i]
        tips.i <- findOS(ancestor = node.i, tree = tree,
                         only.Tip = TRUE, self.include = TRUE,
                         return = "label")
        cNode[i, ] <- apply(data[tips.i, ], 2, fun)
    }
    colnames(cNode) <- colnames(data)
    return(rbind(data, cNode))
}


#' Calculate entity value for each internal node
#'
#' \code{nodeValue.B} calculates value, such as, count, for each internal node.
#' @param tree A phylo object
#' @param data A treeSummarizedExperiment
#' @param fun A function
#'
#' @importFrom utils head
#' @importFrom S4Vectors metadata DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData
#' @keywords internal
#' @return A treeSummarizedExperiment.
#' @author Ruizhu Huang
#'
#'
nodeValue.B <- function(tree, data, fun = sum) {
    if (!inherits(data, "treeSummarizedExperiment")) {
        stop("\n data should be a treeSummarizedExperiment. \n")
    }

    if (!missing(tree)) {
        stop("\n Conflict in tree: phylo object has been found both in data
             and tree argument. \n")
    }
    # extract table and tree from treeSummarizedExperiment.
    tree <- metadata(data)$tree
    table <- assays(data)
    rData <- rowData(data)

    # leaves and internal nodes
    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)

    ## calculate counts for internal nodes
    nN <- length(nodeI)
    nNam <- transNode(tree = tree, input = nodeI)

    # create empty table for values of internal nodes
    tableI <- lapply(table, FUN = function(x){
        y <- matrix(NA, nrow = nN, ncol = ncol(x))
        rownames(y) <- nNam
        colnames(y) <- colnames(data)
        return(y)
    })
    # calculate counts at nodes
    for (i in seq_len(nN)) {
        node.i <- nodeI[i]
        tips.i <- findOS(ancestor = node.i, tree = tree,
                         only.Tip = TRUE, self.include = TRUE,
                         return = "label")
        row.i <- rData$rowID[rData$nodeLab %in% tips.i]
        # calculate count
        for(j in seq_along(tableI)){
            table.j <- table[[j]]
            tableI[[j]][i, ] <- apply(table.j[row.i, ], 2, fun)
        }
    }

    # create empty dataFrame for the rowData of internal nodes
    rDataI <- DataFrame(nodeLab = transNode(tree = tree, input = nodeI,
                                            use.original = TRUE),
                        nodeNum = nodeI,
                        isLeaf = FALSE,
                        rowID = nrow(rData) + seq_len(nN),
                        row.names = nNam )

    # combine rows for internal nodes and leaf nodes
    table <- lapply(table, FUN = function(x){
        rownames(x) <- NULL
        return(x)
    })
    tableI <- lapply(tableI, FUN = function(x){
        rownames(x) <- NULL
        return(x)
    })
    tableC <- Map(rbind, table, tableI)
    rdataC <- rbind(rData, rDataI)

    tse <- treeSummarizedExperiment(tree = tree, assays = tableC,
                                    rowData = rdataC)

    SummarizedExperiment(assays = tableC,
                         rowData = rdataC)
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
#' @param tree A phylo object to provide information of tree structure. We don't
#' need this if data is a treeSummarizedExperiment object.
#' @param data A data frame, matrix or treeSummarizedExperiment object. If it's
#' treeSummarizedExperiment, then tree information is also given.
#' @param fun A function to create the value of an internal node based on values
#' at its descendant leaf nodes. The default is sum.
#'
#' @importFrom utils head
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
#' if(require(ggtree)){
#' (p <- ggtree(tinyTree) + geom_text(aes(label = label)))
#'
#' count <- matrix(rpois(100, 10), nrow =10)
#' rownames(count) <- tinyTree$tip.label
#' colnames(count) <- paste(c("C1_", "C2_"),
#' c(1:5, 1:5), sep = "")
#'
#' count_tinyTree <- nodeValue(data = count,
#' tree = tinyTree, fun = mean)
#'
#' # check to see whether the count of an internal node is the sum
#' # of counts of its descendant leaves.
#' # here, check the first sample as an example
#'
#' nod <- transNode(tree = tinyTree, input = rownames(count_tinyTree))
#' d <- cbind.data.frame(node = nod, count = count_tinyTree[, 1])
#'
#' ggtree(tinyTree) %<+% d + geom_text2(aes(label = count))
#'}
#'
setGeneric("nodeValue", function(tree, data, fun) {
    standardGeneric("nodeValue")
})

#' @rdname nodeValue
#' @export
setMethod("nodeValue", signature("phylo", "matrixDataframe", "function"),
          nodeValue.A)

#' @rdname nodeValue
#' @export
setMethod("nodeValue", signature("ANY", "treeSummarizedExperiment",
                                 "function"), nodeValue.B)

