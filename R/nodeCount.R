#' Calculate the counts at internal nodes
#'
#' \code{nodeCount} is to calculate the counts at internal nodes. The count of an internal node is the sum of counts at its descendant leaves.
#'
#' @param tipTable the tip count table
#' @param wtree a tree (phylo class)
#' @param stree a list of phylo class;
#'              a list of subtrees cut at internal nodes of the wtree
#'
#' @return a count table (matrix class) with a row representing a node and
#' a column representing a sample.
#'
#' @export
#'
#'



nodeCount <- function(tipTable, wtree, stree = NULL) {
    
    if (is.null(tipTable)) {
        stop("tipTable is missing")
    }
    if (!inherits(wtree, "phylo")) {
        stop("object tree is not of class phylo")
    }
    
    # if stree is not provided, generate it using pruneTree
    if (is.null(stree)) {
        cat("stree is not provided and will be generated automatically")
        stree <- pruneTree(tree = wtree)
    } else {
        stree <- stree
    }
    
    # check whether each row of a count table is a tip in the tree
    if (!all(rownames(tipTable) %in% wtree$tip.label)) {
        chx <- sum(!rownames(tipTable) %in% wtree$tip.label)
        warning(chx, " rows could not be found on the tree; only those rows matching the tree tip labels are used")
        tipTable <- tipTable[rownames(tipTable) %in% wtree$tip.label, ]
    }
    
    ## calculate counts for nodes
    nN <- wtree$Nnode
    nNam <- names(stree)
    
    # calculate counts at nodes
    cNode <- matrix(NA, nrow = nN, ncol = ncol(tipTable))
    rownames(cNode) <- nNam
    for (i in 1:nN) {
        node.i <- nNam[i]
        tips.i <- (stree[[node.i]])$tip.label
        cNode[i, ] <- apply(tipTable[tips.i, ], 2, sum)
    }
    colnames(cNode) <- colnames(tipTable)
    return(rbind(tipTable, cNode))
}
