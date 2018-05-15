#' tree aggregation
#'
#' \code{treeAGG} applies minimum p-value algorithm to do tree aggregation
#'
#' @param wtree  a phylo object;
#' @param data data frame (include at least :
#'        1. the label of tree nodes and tips as row names
#'        2. a column of p value
#'        3. a column of adjusted p value
#' @param stree a list of phylo object. The subtrees of \strong{wtree}.
#' @param P.lim the threshold value (for the adjusted p values) to reject a null hypothesis. By default, NULL. If NULL, the algorithm only compares the value provided by \strong{varAGG} and doesn't decide whether to reject a null hypothesis.
#' @param varSIG the name of column for testing significance
#' @param varAGG the name of column used to do tree aggregation
#' @return a vector of character shows the labels of tips and nodes which
#' have selected by the minimum P-value algorithm
#'
#' @export
#'
#'



treeAGG <- function(wtree, data, stree = NULL,
                    P.lim = NULL, varSIG = NULL,
                    varAGG) {

    if (!inherits(wtree, "phylo")) {
        stop("object tree is not of class phylo.")
    }

    # if stree is not provided, generate it using pruneTree
    if (is.null(stree)) {
        cat("stree is not provided and will be generated automatically")
        stree <- pruneTree(tree = wtree)
    } else {
        stree <- stree
    }

    isPhy <- unlist(lapply(stree, FUN = function(x) {
        inherits(x, "phylo")
    }))
    if (!all(isPhy)) {
        stop("object stree is not a list of phylo objects.")
    }

    if (!is.null(P.lim)) {
        if (is.null(varSIG)) {
            stop("varSIG should be specified if P.lim is not null")
        }
    }

    if (!is.null(varSIG)) {
        if (is.null(P.lim)) {
            stop("Specify a value for P.lim if varSIG is not null")
        }
    }
    # tips & nodes from the whole tree
    all.node <- names(stree)
    num.node <- wtree$Nnode
    all.tip <- wtree$tip.label
    num.tip <- length(all.tip)

    ## ----------- tree aggregation----------- firstly, set all tips and nodes as TRUE
    keep.tip <- rep(TRUE, num.tip)
    keep.node <- rep(TRUE, num.node)
    names(keep.tip) <- all.tip
    names(keep.node) <- all.node
    keep <- c(keep.tip, keep.node)

    # compare between the parent and its children.  If the parent has smaller value,
    # set their children as FALSE otherwise set the parent as FALSE.

    for (i in 1:length(all.node)) {
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
            # if min value occur at nodes, remove their children otherwise remove the nodes and
            # keep their children
            if (mRank.node > mRank.child) {
                keep[node.i] <- FALSE
            } else {
                keep[child.i] <- FALSE
            }
        }
    }

    keep1 <- keep[keep]
    namK <- names(keep1)

    if (is.null(P.lim)) {
        final <- namK
    } else {
        isSig <- rownames(data)[data[, varSIG] <= P.lim]
        final <- intersect(isSig, namK)
    }

    return(final)
}
