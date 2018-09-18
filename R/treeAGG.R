
# ----------------------------------------------------------------------------
# if data provided is a data.frame
treeAGG.A <- function(data, sigf.by,
                      sigf.limit, agg.by,
                      tree,
                      message) {

    if (!inherits(tree, "phylo")) {
        stop("object tree is not of class phylo. \n")
    }

    # nodes
    nodeL <- rownames(data)
    nodeN <- transNode(tree = tree, input = nodeL)

    # internal nodes
    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)

    # internal nodes existing in data
    nodeIE <- intersect(nodeN, nodeI)

    # assign keep = TRUE to all nodes
    aggData <- data.frame(nodeNum = nodeN, keep = TRUE)

    # compare between an internal node and its descendant nodes.
    # If the internal node has smaller value, set its descendants as FALSE;
    # otherwise set the internal node as FALSE.

    for (i in seq_along(nodeIE)) {
        node.i <- nodeIE[i]
        desc.i <- findOS(tree = tree, ancestor = node.i,
                         only.Tip = FALSE, self.include = FALSE,
                         return = "number" )
        row.n <- match(node.i, nodeN)
        row.d <- match(desc.i, nodeN)

        # extract the values on the parent and on the children
        agg.n <- data[row.n, agg.by]
        agg.d <- data[row.d, agg.by]

        # replace NA value to 1
        agg.n[is.na(agg.n)] <- 1
        agg.d[is.na(agg.d)] <- 1

        # set FALSE to parent, if the min value is on the children
        # otherwise set FALSE to children
        isMin <- agg.n <= min(agg.d)
        if (isMin) {
            aggData[row.d, "keep"] <- FALSE
        } else {
            aggData[row.n, "keep"] <- FALSE
        }

        if (message) {
            #Sys.sleep(0.001)
            message(i, " out of ", length(nodeIE),
                    " finished", "\r", appendLF = FALSE)
            flush.console()
        }
    }

    final <- data
    final$aggKeep <- aggData$keep

    return(final)
}


# ----------------------------------------------------------------------------
# if data provided is a treeSummarizedExperiment (tse)
treeAGG.B <- function(data, sigf.by,
                      sigf.limit,
                      agg.by, message) {


    # extract tree
    tree <- treeData(data)
    linkD <- linkData(data)

    # leaves and internal nodes
    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)

    # extract row data
    rowD <- rowData(data, internal = FALSE) # orignial rowData
    rowDI <- rowData(data, internal = TRUE) # orignial rowData + result
    namI <- colnames(rowDI)
    indI <- grepl(pattern = "result_assay", x = namI) # result
    rowI <- rowDI[, indI, drop = FALSE]

    # convert to a list
    #  - each element is a data frame whose column is also a data frame
    #  - each element is the result for an element of assays
    rowList1 <- lapply(seq_len(sum(indI)), FUN = function(x) {
        xl <- rowI[, x]
        return(xl)
    })
    names(rowList1) <- colnames(rowI)

    #  - each element is a data frame
    #    it is the result under a specified contrast for an element of assays
    rowList2 <- lapply(rowList1, FUN = function(x) {
        xv <- vector("list", ncol(x))
        names(xv) <- colnames(x)
        for (i in seq_len(ncol(x))) {
            xv[[i]] <- x[, i]
        }
        return(xv)
    })

    # add a column (aggKeep) in each dataframe
    rowList3 <- lapply(rowList2, FUN = function(x) {
        xx <- lapply(x, FUN = function(y) {
            y$aggKeep <- TRUE

            # change the values for nodes with NA value in the result as FALSE
            y$aggKeep[is.na(y[[sigf.by]])] <- FALSE

            return(y)
        })
        return(xx)
    })

    # do tree aggregation by comparing the value from agg.by
    # If the parent node has the smallest value, set the descendants
    # as FALSE; otherwise set the parent node as FALSE.
    for (k  in seq_along(nodeI)) {
        parent.k <- nodeI[k]
        desc.k <- findOS(tree = tree, ancestor = parent.k,
                         only.Tip = FALSE,
                         self.include = FALSE, return = "number")

        # row index
        row.p <- match(parent.k, linkD$nodeNum)
        row.d <- match(desc.k, linkD$nodeNum)

        for (i in seq_along(rowList3)) {
            xi <- rowList3[[i]]
            for (j in seq_along(xi)) {
                xj <- xi[[j]]

                # decide whether to set FALSE to children
                vparent <- xj[[agg.by]][row.p]
                vdesc <- xj[[agg.by]][row.d]

                # input missing values with 1, the maximum value, so that these
                # nodes would not be selected.
                vparent[is.na(vparent)] <- 1
                vdesc[is.na(vdesc)] <- 1

                # If the parent node has the smallest value, set the descendants
                # as FALSE; otherwise set the parent node as FALSE.
                keepParent <- vparent <= min(vdesc)
                if (keepParent) {
                    rowList3[[i]][[j]]$aggKeep[row.d] <- FALSE
                    #xj$aggKeep[row.d] <- FALSE
                } else {
                    rowList3[[i]][[j]]$aggKeep[row.p] <- FALSE
                    #xj$aggKeep[row.p] <- FALSE
                }

                # set FALSE if the sigf.by (e.g. adjusted p-value) value is
                # above sigf.limit (a threshold value) or is NA
                vnode <- xj[[sigf.by]][c(row.p, row.d)]
                notSigf <- (vnode > sigf.limit | is.na(vnode))
                row.ns <- c(row.p, row.d)[notSigf]
                rowList3[[i]][[j]]$aggKeep[row.ns] <- FALSE
            }
        }

        if (message) {
            message(k, " out of ", length(nodeI),
                    " is done. ", "\r", appendLF = FALSE)
            flush.console()
        }
    }

    # reshape: convert a list into a dataFrame
    rowList4 <- lapply(seq_along(rowList3), FUN = function(j) {
            x <- rowList3[[j]]
            dx <- x[[1]][, 0]
            for (i in seq_along(x)) {
                xi <- x[[i]]
                nxi <- names(x)[i]
                dx[[nxi]] <- xi

            }
            return(dx)
        })
    names(rowList4) <- names(rowList3)

    # update rowData
    res <- data
    rowData(res) <- rowData(data, internal = FALSE)

    for (i in seq_along(rowList4)) {
       rowData(res)[[names(rowList4)[i]]] <- as(rowList4[[i]], "internal_rowData")
    }

   return(res)

}


# ----------------------------------------------------------------------------
#' Tree aggregation
#'
#' \code{treeAGG} combines the p values with the tree structure and decide the
#' which nodes to be aggregated to based on the min-p algorithm.
#'
#' @param data A data frame or a treeSummarizedExperiment.
#'      \itemize{
#'      If a data frame, it should include at least:
#'      \item a column of node labels
#'           (use labels from this column to map each row to a node of tree.)
#'      \item a column for tree aggregation
#'           (use value from this column to decide whether to aggregate.)
#'      \item a column of adjusted p value
#'           (use value from this column to decide whether to reject a null
#'           hypothesis.)
#'      }
#' @param sigf.by A column name. The column contains the p value or adjusted
#' p value.
#' @param agg.by A column name. The column used to do tree aggregation.
#' Commonly, it is the column including p value or adjusted p value.
#' @param sigf.limit A numeric value. The threshold value (for p value or
#'   adjusted p value) to reject a null hypothesis. The chosen value depends on
#'   the \code{sigf.by}.
#' @param tree  A phylo object. A optional argument. Only use when \code{data}
#'   is a data frame.
#' @param message A logical value. The default is TRUE. If TRUE, it will print
#'   out the currenet status of a process.
#' @importFrom S4Vectors DataFrame
#' @importFrom utils flush.console
#'
#' @return A data frame
#' @author Ruizhu Huang
#' @name treeAGG
#' @export
#' @examples
#'
#' # We recommend to use treeAGG as example 2. It works in the case of example 1.
#' # However, it is difficult to rownames a data frame when some internal nodes
#' # of the tree doesn't have labels. Furthermore, it could do tree agg
#' # Example 1
#' library(ggtree)
#' data(tinyTree)
#'
#' # data
#' set.seed(3)
#' pv <- runif(19)
#' pValue <- rank(pv)/length(pv)*pv
#' treeLab <- c(tinyTree$tip.label, tinyTree$node.label)
#' df <- cbind.data.frame(pV = pValue,
#' stringsAsFactors = FALSE)
#' rownames(df) <- treeLab
#'
#'
#'
#' # tree aggregation
#' (tt <- treeAGG(tree = tinyTree, data = df, sigf.limit = 0.05,
#' sigf.by = "pV", agg.by = "pV"))
#'
#' # display the tree structure and p value at each node
#' tt$node <- transNode(tree = tinyTree, input = rownames(tt))
#'
#' # p value at each node is given as blue number in tree
#' # the selected nodes after aggregation is labelled with orange points
#' # these selected nodes have lower p-value than its descendant nodes if they
#' # have descendant nodes.
#' ggtree(tinyTree) %<+% tt + geom_text2(aes(label = label), hjust = -0.2) +
#' geom_text2(aes(label = round(pv, 3)), vjust = -0.5, color = "blue",
#'  hjust = -0.15) +
#' geom_point2(aes(subset = aggKeep), color = "orange", size = 2)
#'
#'
#' # Example 2
#' set.seed(1)
#' y <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(y) <- paste(rep(LETTERS[1:3], each = 10), rep(1:10,3), sep = "_")
#' rownames(y) <- tinyTree$tip.label
#'
#' rowInf <- data.frame(nodeLab = rownames(y),
#'                     var1 = sample(letters[1:3], 10, replace = TRUE),
#'                     var2 = sample(c(TRUE, FALSE), 10, replace = TRUE),
#'                     stringsAsFactors = FALSE)
#' colInf <- data.frame(gg = factor(sample(1:3, 30, replace = TRUE)),
#'                     group = rep(LETTERS[1:3], each = 10))
#' toy_lse <- leafSummarizedExperiment(tree = tinyTree,
#'                                     assays = list(y, (2*y), 3*y),
#'                                     rowData = rowInf,
#'                                     colData = colInf)
#'
#' toy_tse <- nodeValue(data = toy_lse, fun = sum, message = TRUE)
#'
#' new_tse <- runEdgeR(obj = toy_tse, use.assays = 1, design = NULL,
#'                     contrast = NULL, normalize = TRUE, method = "TMM",
#'                     adjust.method = "BH")
#' # option 1: provide treeSummarizedExperiment
#' outR1 <- treeAGG(data = new_tse)
#'
#'  # or we could extract the result and the tree structure to do aggregation
#'  (res <- rowData(new_tse)$result_assay1$contrastNULL)
#'  rownames(res) <- linkData(new_tse)$nodeLab
#'  (Tree <- treeData(new_tse))
#'
#'  # option2: provide data frame and tree structure
#' outR2 <- treeAGG(data = res, tree = Tree)
#'
#'
#'
#'
setGeneric("treeAGG", function(data, sigf.by = "FDR",
                               sigf.limit = 0.05, agg.by = "FDR",
                               tree, message  = FALSE) {
    standardGeneric("treeAGG")
})

#' @rdname treeAGG
#' @importFrom S4Vectors DataFrame
setMethod("treeAGG", signature(data = "treeSummarizedExperiment"),
          treeAGG.B)

#' @rdname treeAGG
setMethod("treeAGG", signature(data = "data.frame"),
          treeAGG.A)

#' @rdname treeAGG
#' @importClassesFrom S4Vectors DataFrame
setMethod("treeAGG", signature(data = "DataFrame"),
          treeAGG.A)








