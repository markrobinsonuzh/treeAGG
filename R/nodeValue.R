# when data is a data frame or a matrix
nodeValue.A <- function(data, fun = sum,
                        tree, message = FALSE,
                        level = NULL) {
    if (!(is.data.frame(data) |
          is.matrix(data) | is(data, "DataFrame")) ) {
        stop("data should be a matrix or data.frame")
    }

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (is.null(data[["nodeLab_alias"]])) {
        lab1 <- data[["nodeLab"]]
        if (is.null(lab1)) {
            lab1 <- rownames(data)
        }
    }
    num1 <- transNode(tree = tree, input = lab1)
    num2 <- unique(as.vector(tree$edge))
    if (!all(num1 %in% num2)) {
        chx <- setdiff(num1, num2)
        chx <- transNode(tree = tree, input = chx)
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
        nodeA <- level
        }

    ## Find the rows of descendants
    nodeA.l <- intersect(nodeA, leaf)
    nodeA.i <- intersect(nodeA, nodeI)

    if (length(nodeA.i)) {
        desI <- findOS(tree = tree, ancestor = nodeA.i, only.leaf = TRUE,
                       self.include = FALSE, use.alias = TRUE,
                       message = message)
    } else {
        desI <- numeric(0)
    }

    desL <- as.list(nodeA.l)
    desA <- c(desL, desI)
    leafRow <- lapply(desA, FUN = function(x){
        # find the rows of descendant leaves in the assay table
        # Note: don't use match!!! (duplicated values might exist in tipNum)
        ri <- which(tipNum %in% x)
        return(ri)
    })

    ## generate values for nodes from their descendant leaves
    # assays
    ## Extract the colum used to do aggregation from the link data

    listD1 <- lapply(leafRow, FUN = function(x){
        tabA[x, , drop = FALSE]
    })

    listD2 <- lapply(listD1, FUN = function(x){
        as.list(x)
    })

    listD3 <- rapply(listD2, f = fun, how = "list")

    listD4 <- lapply(listD3, FUN = function(x) {
        xx <- do.call(cbind.data.frame, c(x, stringsAsFactors = FALSE))
        colnames(xx) <- names(x)
        xx
    })

    tabN <- do.call(rbind, listD4)

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
    if (ncol(rData) == 0) {
        rData <- DataFrame(nodeLab = rownames(data))
    }
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


    ## Find the rows of descendants
    leafRow <- lapply(seq_along(nodeA), FUN = function(x) {
        # get a node
        xx <- nodeA[x]

        # find its descendant leaves
        lx <- findOS(tree = tree, ancestor = xx,
                     only.leaf = TRUE, self.include = TRUE,
                     use.alias = TRUE)

        # find the rows of descendant leaves in the assay table
        # Note: don't use match!!! (duplicated values might exist in tipNum)
        ri <- which(tipNum %in% lx)
        rx <- tipNum[ri]
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
    rdataA <- do.call(rbind.data.frame, rD)
    rownames(rdataA) <- NULL
    colnames(rdataA) <- colnames(rData)
    rdataA <- DataFrame(rdataA)
    rdataA$nodeLab <- transNode(tree = tree, input = nodeA,
                                use.alias = TRUE, message = FALSE)


    mData.new <- mData[!names(mData) %in% "tree"]
    tse <- treeSummarizedExperiment(tree = tree,
                                    assays = tabNA, metadata = mData.new,
                                    colData = cData, rowData = rdataA)
    return(tse)


}

# when data is a treeSummarizedExperiment
nodeValue.C <- function(data, fun = sum, message = FALSE,
                        level = NULL, col.linkData = "nodeNum") {

    if (!is(data, "treeSummarizedExperiment")) {
        stop("\n data should be a leafSummarizedExperiment. \n")
    }

    # =========================================================================
    if (message) {cat("Determine the level on the tree...")}
    ## Extract the colum used to do aggregation from the link data
    lk <- linkData(data)
    av <- lk[[col.linkData]]

    ## Decide the knots
    if (!is.null(level)) {
        isOut <- !level %in% av
        if (any(isOut)) {
            stop(level[isOut], " can't be found in the specified column ",
                 col.linkData, "\n")
        }

        # the value indexed
        ind <- which(av %in% level)
    } else {
        ind <- seq_along(av)
     }

    av <- av[ind]
    nd <- lk[["nodeNum"]][ind]


    if (col.linkData %in% c("nodeNum", "nodeLab", "nodeLab_alias")) {
        lnd <- findOS(tree = treeData(data), ancestor = nd,
                          only.leaf = TRUE, self.include = TRUE,
                          use.alias = TRUE, message = FALSE)
        knot <- as.list(level)
    } else {
        lnd <- base::split(x = nd, f = av)
        knot <- lapply(seq_along(lnd), FUN = function(x) {
            if (message) {cat(x, " out of ", length(lnd), "\n")}
            vx <- lnd[[x]]
            if (length(vx)) {
                xx <- signalNode(tree = treeData(data), node = vx,
                                 use.alias = FALSE)
            } else {
                xx <- numeric(0)
            }
            return(xx)
        })
    }

    idK <- lapply(seq_along(lnd), FUN = function(x) {
        ix <- lnd[[x]]
        ik <- lk[["rowID"]][ix]
        unlist(ik)
    })
    # =========================================================================
    if (message) { cat("Aggregate the tables in assays... \n")}
    ## Do aggregation on tables of assays
    # assays
    tabA <- assays(data, use.nodeLab = TRUE)
    aggK <- vector("list", length = length(tabA))
    for (i in seq_along(tabA)) {
        tabL.i <- lapply(idK, FUN = function(x) {
            tabA.i <- tabA[[i]]
            xi <- tabA.i[x, , drop = FALSE]
            ri <- apply(xi, 2, fun)
            return(ri)
        })
        tab.i <- do.call(rbind, tabL.i)
        aggK[[i]] <- tab.i

    }

    # =========================================================================
    if (message) { cat("Aggregate the row data... \n")}
    ## do aggregation on rowData
    rData <- rowData(data)
    if (ncol(rData)) {
        listD1 <- lapply(idK, FUN = function(x){
            rData[x, , drop = FALSE]
        })

        listD2 <- lapply(listD1, FUN = function(x){
            as.list(x)
        })

        listD3 <- rapply(listD2, f = function(x) {
            xx <- unique(x)
            ifelse(length(xx) > 1, NA, xx)
        }, how = "list")

        listD4 <- lapply(listD3, FUN = function(x) {
            xx <- do.call(cbind.data.frame, x)
            colnames(xx) <- names(x)
            xx
        })

        rDataA <- do.call(rbind, listD4)
        rownames(rDataA) <- NULL
    } else {
        rDataA <- rData[rep(1, nrow(aggK[[1]])), ]
    }
    rDataA <- DataFrame(rDataA)
    rownames(rDataA) <- names(lnd)
    # =========================================================================
    if (message) { cat("Aggregate the link data ... \n")}
    ## do aggregation on linkData
    # lkData <- linkData(data)
    # lkData <- lkData[lkData$isLeaf, ]
    # if (is.null(lkData$nodeLab_alias)) {
    #     rownames(lkData) <- lkData$nodeLab
    # } else {
    #     rownames(lkData) <- lkData$nodeLab_alias
    # }
    #
    # treeD <- treeData(data)
    # nkData <- nodeValue.A(data = data.frame(lkData),
    #                       fun = function(x) {
    #     ux <- unique(x)
    #     ux <- as.charactor(ux)
    #     ifelse(length(ux) == 1, ux, NA)
    # }, tree = treeD, message = message, level = nd)



    if (message) { cat("Output the new treeSummarizedExperiment... \n")}

    lenK <- unlist(lapply(knot, length))
    if (any(lenK > 1)) {
        # Give warning if a knot includes more than one branch and output
        # results as a summarizedExperiment object
        ik <- lenK > 1
        iav <- unique(av[which(ik %in% TRUE)])
        warning(paste(iav, collapse = " "),
                ": each has nodes on different branches. \n")
        out <- SummarizedExperiment(assays = aggK, metadata = metadata(data),
                                    colData = colData(data), rowData = rDataA)
        warning("The output is a SummarizedExperiment object.")
    } else {
        rDataA$nodeLab <- transNode(tree = treeData(data), input = unlist(knot),
                                    use.alias = TRUE, message = FALSE)

        out <- treeSummarizedExperiment(tree = treeData(data),
                                        linkData = linkData(data),
                                        assays = aggK, metadata = metadata(data),
                                        colData = colData(data), rowData = rDataA)
    }

    return(out)

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
#' @param ... The additional argument allows only for
#'   \code{treeSummarizedExperiment}.
#' @param col.linkData the column names from \code{linkData}.
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
                                 level = NULL, ...) {
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

