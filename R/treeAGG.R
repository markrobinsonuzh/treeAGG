
# ----------------------------------------------------------------------------
# if data provided is a data.frame
treeAGG.A <- function(data, sigf.by = NULL,
                      sigf.limit = 0.05, agg.by = NULL,
                      tree = NULL, node.by = "nodeLab") {

    if (!inherits(tree, "phylo")) {
        stop("object tree is not of class phylo. \n")
    }

    if (!inherits(data[, node.by], "character")) {
        stop("The column ", node.by, " should be character. \n")
    }
    # nodes
    nodeL <- data[, node.by]
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
    }

    final <- data
    final$aggKeep <- aggData$keep

    return(final)
}


# ----------------------------------------------------------------------------
# if data provided is a treeSummarizedExperiment (tse)
treeAGG.B <- function(data, sigf.by = NULL,
                        sigf.limit = 0.05, agg.by = NULL) {


    # extract tree
    tree <- treeData(data)
    linkD <- linkData(data)
    assayD <- assays(data)

    # leaves and internal nodes
    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)


    # assign keep = TRUE to all nodes
    keep <- DataFrame(nodeLab = linkD$nodeLab,
                      nodeNum = linkD$nodeNum, aggKeep = TRUE)
    keepList <- lapply(seq_along(assayD), FUN = function(x){
        cbind(assayD[[x]], keep)
    })

    # extract a column for aggregation and a column for rejecting hypothesis
    # test
    aggList<- lapply(assayD, FUN = function(x){
        data.frame(nodeNum = linkD$nodeNum,
                   varSig = x[linkD$rowID, sigf.by],
                   varAgg = x[linkD$rowID, agg.by])
    })

    # compare between an internal node and its descendant nodes.
    # If the internal node has smaller value, set its descendants as FALSE;
    # otherwise set the internal node as FALSE.
    nodeIE <- intersect(nodeI, linkD$nodeNum)
    for (i in seq_along(nodeIE)) {
        # an node (parent) & its descendants (children)
        node.i <- nodeIE[i]
        desc.i <- findOS(tree = tree, ancestor = node.i, only.Tip = FALSE,
                         self.include = FALSE, return = "number")
        # row index
        row.n <- match(node.i, linkD$nodeNum)
        row.d <- match(desc.i, linkD$nodeNum)

        # decide whether to set FALSE to children
        dc.i <- lapply(seq_along(keepList), FUN = function(x) {
            vnode <- aggList[[x]]$varAgg[row.n]
            vdesc <- aggList[[x]]$varAgg[row.d]

            # set NA value to 1 (to avoid the NA p-value caused by filteration
            # or not observed)
            vnode[is.na(vnode)] <- 1
            vdesc[is.na(vdesc)] <- 1

            # the min value at parent node (node.i)
            vnode <= min(vdesc)

        })

        # if decision (stored in dc.i) is true, set FALSE to children
        # otherwise, set FALSE to the parent.
        for (j in seq_along(keepList)) {
            if (dc.i[[j]]) {
                keepList[[j]]$aggKeep[row.d] <- FALSE
            } else {
                keepList[[j]]$aggKeep[row.n] <- FALSE
            }
        }

    }

    # set FALSE if the varSig (e.g. adjusted p-value) value is above sigf.limit
    # (a threshold value)
    nodeA <- linkD$nodeNum
    for (j in seq_along(keepList)) {
        # find the pair of rows that have the same nodeNum in aggList
        # and keepList
        sigJ <- aggList[[j]]$varSig
        nodeJ <- aggList[[j]]$nodeNum
        for (i in seq_along(nodeA)) {
            node.i <- nodeA[i]
            ind.i <- nodeJ == node.i
            sig.i <- sigJ[ind.i]
            if (sig.i > sigf.limit) {
                keepList[[j]]$aggKeep[i] <- FALSE
            }
        }
    }

    return(keepList)


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
#' @param node.by A column name. The column stores the node label. A optional
#'   argument. Only use when \code{data} is a data frame.


#'
#' @return A data frame
#' @author Ruizhu Huang
#' @name treeAGG
#' @export
#' @examples
#'
#' library(ggtree)
#'
#' data(tinyTree)
#'
#' # data
#' set.seed(3)
#' pv <- runif(19)
#' pValue <- rank(pv)/length(pv)*pv
#' treeLab <- c(tinyTree$tip.label, tinyTree$node.label)
#' df <- cbind.data.frame(pV = pValue,
#' stringsAsFactors = FALSE, nodeLab = treeLab)
#'
#'
#' # tree aggregation
#' (tt <- treeAGG(tree = tinyTree, data = df, sigf.limit = 0.05,
#' sigf.by = "pV", agg.by = "pV", node.by ="nodeLab"))
#'
#' # display the tree structure and p value at each node
#' tt$node <- transNode(tree = tinyTree, input = tt$nodeLab)
#'
#' # p value at each node is given as blue number in tree
#' # the selected nodes after aggregation is labelled with orange points
#' # these selected nodes have lower p-value than its descendant nodes if they
#' # have descendant nodes.
#' ggtree(tinyTree) %<+% tt + geom_text2(aes(label = label), hjust = -0.2) +
#' geom_text2(aes(label = round(pV, 3)), vjust = -0.5, color = "blue",
#'  hjust = -0.15) +
#' geom_point2(aes(subset = aggKeep), color = "orange", size = 2)
#'
#'
#'
setGeneric("treeAGG", function(data, sigf.by,
                               sigf.limit, agg.by,
                               tree, node.by) {
    standardGeneric("treeAGG")
})

#' @rdname treeAGG
setMethod("treeAGG", signature(data = "treeSummarizedExperiment"),
          treeAGG.B)

#' @rdname treeAGG
setMethod("treeAGG", signature(data = "data.frame"),
          treeAGG.A)

#' @rdname treeAGG
#' @importClassesFrom S4Vectors DataFrame
setMethod("treeAGG", signature(data = "DataFrame"),
          treeAGG.A)








