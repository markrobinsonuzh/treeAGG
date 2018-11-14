# when data is a data frame or a matrix
nodeValue.A <- function(data, fun = sum,
                        tree, message = FALSE,
                        level = NULL) {
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
    ## rename data with the alias of the node labels
    tabA <- data
    rn <- rownames(tabA)
    tipNum <- transNode(tree = tree, input = rn, message = FALSE)
    leafLab_alias <- transNode(tree = tree, input = tipNum,
                               use.alias = TRUE, message = FALSE)
    rownames(tabA) <- leafLab_alias

    ## nodes
    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)
    if (is.null(level)) {
        nodeA <- c(leaf, nodeI)
    } else {
        if (is.character(level)) {
            level <- transNode(tree = tree, input = level,
                               use.alias = FALSE, message = FALSE)
        }
        ll <- level
        lo <- .findParallel(tree = tree, input = level, use.alias = FALSE)
        nodeA <- c(ll, lo)
    }

    ## Find the rows of descendants
    leafRow <- lapply(seq_along(nodeA), FUN = function(x) {
        # get a node
        xx <- nodeA[x]

        # find its descendant leaves
        lx <- findOS(tree = tree, ancestor = xx,
                     only.Tip = TRUE, self.include = TRUE,
                     use.alias = TRUE)

        # find the rows of descendant leaves in the assay table
        rx <- match(lx, tipNum)
        return(rx)
    })

    ## generate values for nodes from their descendant leaves
    # assays
    tabNA <- lapply(leafRow, FUN = function(x) {
        tab.x <- tabA[x, , drop = FALSE]
        rx <- apply(tab.x, 2, fun)
        rx
    })

    tabN <- do.call(rbind, tabNA)

    ## the node labels
    nodeLab <- transNode(tree = tree, input = nodeA,
                         use.alias = FALSE, message = FALSE)
    if (anyDuplicated(nodeLab)) {
        nodeLab_alias <- transNode(tree = tree, input = nodeA,
                                   use.alias = TRUE, message = FALSE)
        rownames(tabN) <- nodeLab_alias
    } else {
        rownames(tabN) <- nodeLab
    }


    # output
    return(tabN)
}

# when data is a leafSummarizedExperiment
nodeValue.B <- function(data, fun = sum, message = FALSE,
                        level = NULL) {
    if (!is(data, "leafSummarizedExperiment")) {
        stop("\n data should be a leafSummarizedExperiment. \n")
    }


    # extract table and tree from treeSummarizedExperiment.
    tree <- metadata(data)$tree
    tabA <- assays(data)
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
    if (is.null(level)) {
        nodeA <- c(leaf, nodeI)
    } else {
        if (is.character(level)) {
            level <- transNode(tree = tree, input = level,
                               use.alias = FALSE, message = FALSE)
        }
        ll <- level
        lo <- .findParallel(tree = tree, input = level, use.alias = FALSE)
        nodeA <- c(ll, lo)
    }

    nN <- length(nodeA)
    ## Find the rows of descendants
    leafRow <- lapply(seq_along(nodeA), FUN = function(x) {
        # get a node
        xx <- nodeA[x]

        # find its descendant leaves
        lx <- findOS(tree = tree, ancestor = xx,
                     only.Tip = TRUE, self.include = TRUE,
                     use.alias = TRUE)

        # find the rows of descendant leaves in the assay table
        rx <- match(lx, tipNum)
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

    # Create the rowData
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

    ## Check whether there are duplicated values in nodeLab
    # If yes, add another column nodeLab_alias
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

# when data is a treeSummarizedExperiment
nodeValue.C <- function(data, fun = sum, message = FALSE,
                        level = NULL) {
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

    if (is.null(level)) {
        nodeA <- c(leaf, nodeI)
    } else {
        if (is.character(level)) {
            level <- transNode(tree = tree, input = level,
                               use.alias = FALSE, message = FALSE)
        }
        ll <- level
        lo <- .findParallel(tree = tree, input = level, use.alias = FALSE)
        nodeA <- c(ll, lo)
    }


    ## all nodes
    nN <- length(nodeA)

    # find the rows of descendants
    leafRow <- lapply(seq_along(nodeA), FUN = function(x) {
        # get nodes
        y <- nodeA[x]

        # find the descendant leaves
        lx <- findOS(tree = tree, ancestor = y,
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
#' @param level A vector of nodes. This indicates the level to which the
#'   aggregation is made. For example, a tree has two main branches (A & B).
#'   Each branch (A or B) has multiple sub-branches. If the aggregation to the
#'   two main branches is desired, we could use \code{ level = c(A, B)},
#'   \code{level = A} or \code{level = B}. In other words, users could specify
#'   the level by provide some branches that covers parts of the tree. The rest
#'   part of the tree is pruned automactially with the minimum number of
#'   branches.
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
#' ## The node labels are in orange and the node numbers are in blue
#' (ggtree(tinyTree,branch.length = 'none')+
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7))

#' set.seed(1)
#' count <- matrix(rpois(100, 10), nrow =10)
#' rownames(count) <- tinyTree$tip.label
#' colnames(count) <- paste(c("C1_", "C2_"),
#' c(1:5, 1:5), sep = "")
#'
#' ## if the level isn't specified, the aggregation is made at each level
#' (agg1 <- nodeValue(data = count, tree = tinyTree,
#'                   fun = sum, level = NULL))
#'
#' ## specify the level that the aggregation should be made to
#'  # 13 (a branch); 16 (a branch); 9 & 10 (are automatically found and
#'  #  separate because they are not in a same branch)
#' (agg2 <- nodeValue(data = count, tree = tinyTree,
#'                   fun = sum, level = c(13, 16)))
#'
#'  # 13 (a branch); 15, 10 (are automatically found and
#'  #  separate because they are not in a same branch)
#' (agg3 <- nodeValue(data = count, tree = tinyTree,
#'                   fun = sum, level = c(13)))
#'
#' # check to see whether the count of an internal node is the sum
#' # of counts of its descendant leaves.
#' # here, check the first sample as an example
#'
#' nod <- transNode(tree = tinyTree, input = rownames(agg1),
#'                  message = FALSE)
#' d <- cbind.data.frame(node = nod, count = agg1[, 1])
#'
#' ggtree(tinyTree) %<+% d + geom_text2(aes(label = count))
#'}
#'
setGeneric("nodeValue", function(data, fun = sum,
                                 tree, message = FALSE,
                                 level = NULL) {
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
          signature(data = "leafSummarizedExperiment"),
          nodeValue.B)

#' @rdname nodeValue
#' @importFrom utils flush.console
#' @importFrom methods is
setMethod("nodeValue",
          signature(data = "treeSummarizedExperiment"),
          nodeValue.C)

