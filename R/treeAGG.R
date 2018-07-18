#' Tree aggregation
#'
#' \code{treeAGG} combines the p values with the tree structure and decide the
#' which nodes to be aggregated to based on the min-p algorithm.
#'
#' @param tree  A phylo object;
#' @param data A data frame (include at least :
#'        1. the label of tree nodes and tips as row names
#'        2. a column for tree aggregation
#'           (use value from this column to decide whether to aggregate)
#'        3. a column of adjusted p value
#'           (use value from this column to decide whether to reject a null
#'           hypothesis)
#' @param stree A list of phylo object. The subtrees of \strong{tree}.
#' @param pLim The threshold value (for \strong{varSIG}) to reject a null
#'   hypothesis. By default, NULL. If NULL, the algorithm only compares the
#'   value provided by \strong{varAGG} and doesn't decide whether to reject a
#'   null hypothesis.
#' @param varSIG The name of column storing adjusted p values.
#' @param varAGG The name of column used to do tree aggregation, eg. the name of
#'   the p value or adjusted p value column
#'
#' @export
#' @return A character vector, containing the labels of tips and nodes which
#' have selected by the minimum P-value algorithm
#' @author Ruizhu Huang
#'
#' @examples
#'
#' library(ggtree)
#'
#' data(tinyTree)
#'
#' # data
#' set.seed(3)
#' pv <- runif(19)
#' apv <- rank(pv)/length(pv)*pv
#' df <- cbind.data.frame(pval = pv, adj_p = apv)
#' rownames(df) <- c(tinyTree$tip.label, tinyTree$node.label)
#'
#' # display the tree structure and p value at each node
#' df1 <- cbind.data.frame(
#' node = transNode(tree = tinyTree, input = rownames(df)),
#' df)
#'
#' ggtree(tinyTree) %<+% df1 + geom_text2(aes(label = label)) +
#' geom_text2(aes(label = round(adj_p, 3)), vjust = -0.5, color = "blue")
#'
#' # tree aggregation
#' (tt <- treeAGG(tree = tinyTree, data = df, pLim = 0.05,
#' varSIG = "adj_p", varAGG = "pval"))
#'


treeAGG <- function(tree, data, stree = NULL,
                    pLim = NULL, varSIG = NULL,
                    varAGG) {

    if (!inherits(tree, "phylo")) {
        stop("object tree is not of class phylo.")
    }

    # if stree is not provided, generate it using pruneTree
    if (is.null(stree)) {
        #cat("stree is not provided and will be generated automatically")
        stree <- pruneTree(tree = tree)
    } else {
        stree <- stree
    }

    isPhy <- unlist(lapply(stree, FUN = function(x) {
        inherits(x, "phylo")
    }))
    if (!all(isPhy)) {
        stop("object stree is not a list of phylo objects.")
    }

    if (!is.null(pLim)) {
        if (is.null(varSIG)) {
            stop("varSIG should be specified if pLim is not null")
        }
    }

    if (!is.null(varSIG)) {
        if (is.null(pLim)) {
            stop("Specify a value for pLim if varSIG is not null")
        }
    }
    # tips & nodes from the whole tree
    all.node <- names(stree)
    num.node <- tree$Nnode
    all.tip <- tree$tip.label
    num.tip <- length(all.tip)

    ## ----------- tree aggregation----------- firstly, set all tips and nodes as
    TRUE
    keep.tip <- rep(TRUE, num.tip)
    keep.node <- rep(TRUE, num.node)
    names(keep.tip) <- all.tip
    names(keep.node) <- all.node
    keep <- c(keep.tip, keep.node)

    # compare between the parent and its children.  If the parent has smaller
    # value, set their children as FALSE otherwise set the parent as FALSE.

    for (i in seq_along(all.node)) {
        # parent & children
        node.i <- all.node[i]  # parent node
        tree.i <- stree[[node.i]]
        child.i <- setdiff(c(tree.i$node.label, tree.i$tip.label), node.i)

        # add 1 here to avoid the NA p-value in tips / nodes (NA might be due to the
        # filteration in DESeq or not observed)
        rank <- data[, varAGG]
        names(rank) <- rownames(data)
        mRank.child <- min(c(rank[child.i], 1), na.rm = TRUE)
        mRank.node <- min(c(rank[node.i], 1), na.rm = TRUE)
        if (keep[node.i]) {
            # if min value occur at nodes, remove their children otherwise remove the
            # nodes and keep their children
            if (mRank.node > mRank.child) {
                keep[node.i] <- FALSE
            } else {
                keep[child.i] <- FALSE
            }
        }
    }

    keep1 <- keep[keep]
    namK <- names(keep1)

    if (is.null(pLim)) {
        final <- namK
    } else {
        isSig <- rownames(data)[data[, varSIG] <= pLim]
        final <- intersect(isSig, namK)
    }

    return(final)
    # datF <- cbind.data.frame(label = rownames(data),
    #                         select = ifelse(rownames(data) %in% final,
    #                                         TRUE, FALSE))
    #return(datF)

}
