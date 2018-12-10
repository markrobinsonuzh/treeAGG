
# when data is a data frame or a matrix
nodeValue.A <- function(data, fun = sum, tree, message = FALSE) {
    if (!(is.data.frame(data) |
          is.matrix(data))) {
        stop("data should be a matrix or data.frame")
    }

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (!setequal(rownames(data), tree$tip.label)) {
        chx <- setdiff(rownames(data), tree$tip.label)
        chs <- head(chx)
        stop(cat("Some row names can't be matched to the labels of leaf nodes:",
                 chs, "\n"))
    }
    ## rename data with the alias of the node label
    rn <- rownames(data)
    leafNum <- transNode(tree = tree, input = rn, message = FALSE)
    leafLab_alias <- transNode(tree = tree, input = leafNum,
                               use.alias = TRUE, message = FALSE)
    rownames(data) <- leafLab_alias

    ## nodes
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
                         use.alias = TRUE)
        tips.i <- names(tips.i)
        cNode[i, ] <- apply(data[tips.i, ], 2, fun)

        # print out the running process
        if (message) {
            message(i, " out of ", nN , " finished", "\r", appendLF = FALSE)
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


# when data is a TreeSummarizedExperiment
nodeValue.B <- function(data, fun = sum, message = FALSE) {
    if (!is(data, "treeSummarizedExperiment")) {
        stop("\n data should be a leafSummarizedExperiment. \n")
    }


    # extract table and tree from treeSummarizedExperiment.
    tree <- treeData(data)
    tabA <- assays(data, use.nodeLab = TRUE)
    rData <- rowData(data, use.names = FALSE)
    lData <- linkData(data)
    cData <- colData(data)
    mData <- metadata(data)

    # leaf nodes and internal nodes
    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)
    nodeA <- c(leaf, nodeI)

    ## all nodes
    nN <- length(nodeA)

    # find the rows of descendants
    leafRow <- lapply(seq_along(nodeA), FUN = function(x) {
        # get a node
        xx <- nodeA[x]

        # find its descendant leaves
        lx <- findOS(tree = tree, ancestor = xx,
                     only.Tip = TRUE, self.include = TRUE,
                     use.alias = TRUE)

        # find the rows of descendant leaves in the assay table
        ind <- match(lx, lData[["nodeNum"]])
        rx <- lData[["rowID"]][ind]
        return(rx)
    })

    ## generate values for nodes from their descendant leaves
    # assays
    tabNA <- vector("list", length = length(tabA))

    for (i in seq_along(tabA)) {
        tabL.i <- lapply(leafRow, FUN = function(x) {
            tabA.i <- tabA[[i]]
            xi <- tabA.i[x, , drop = FALSE]
            ri <- apply(xi, 2, fun)
            return(ri)
        })
        tab.i <- do.call(rbind, tabL.i)
        tabNA[[i]] <- tab.i

    }

    # rowData
    rD <- lapply(leafRow, FUN = function(x) {
        xi <- rData[x, , drop = FALSE]
        ri <- apply(xi, 2, FUN = function(x) {
            ui <- unique(x)
            si <- ifelse(length(ui) > 1, NA, ui)
            return(si)
        })
        return(ri)
    })
    rdataA <- do.call(rbind, rD)
    rownames(rdataA) <- NULL

    ## Create the node labels
    nodeLab <- transNode(tree = tree, input = nodeA,
                         use.alias = FALSE, message = FALSE)

    if (anyDuplicated(nodeLab)) {
        nodeLab_alias <- transNode(tree = tree, input = nodeA,
                                   use.alias = TRUE, message = FALSE)
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


    mData.new <- mData
    tse <- treeSummarizedExperiment(linkData = linkD, tree = tree,
                                    assays = tabNA, metadata = mData.new,
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
#' @return If input data is a treeSummarizedExperiment object, a
#'   treeSummarizedExperiment is returned; otherwise, a matrix is returned
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
setGeneric("nodeValue", function(data, fun = sum,
                                 tree, message = FALSE) {
    standardGeneric("nodeValue")
})

#' @rdname nodeValue
#' @importFrom utils flush.console
setMethod("nodeValue", signature(data = "ANY"),
          nodeValue.A)

#' @rdname nodeValue
#' @importFrom utils flush.console
#' @importFrom methods is
setMethod("nodeValue",
          signature(data = "treeSummarizedExperiment"),
          nodeValue.B)
