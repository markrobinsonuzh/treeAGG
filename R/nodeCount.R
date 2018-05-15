#' Calculate the counts at internal nodes
#'
#' \code{nodeCount} is to calculate the counts at internal nodes. The count of an internal node is the sum of counts at its descendant leaves.
#'
#' @param tipTable the tip count table
#' @param wtree a tree (phylo class)
#' @param stree a list of phylo class;
#'              a list of subtrees cut at internal nodes of the wtree
#' @param fun a function to create the count of an internal node based on the counts at its descendant leaf nodes. The default is sum
#' @return a count table (matrix class) with a row representing a node and
#' a column representing a sample.
#'
#' @export
#'
#' @examples
#' data("tinyTree")
#' (p <- ggtree(tinyTree) + geom_text(aes(label = label)))
#'
#' count <- matrix(rpois(100, 10), nrow =10)
#' rownames(count) <- tinyTree$tip.label
#' colnames(count) <- paste(c("C1_", "C2_"),
#' c(1:5, 1:5), sep = "")
#'
#' count_tinyTree <- nodeCount(tipTable = count,
#' wtree = tinyTree, fun = mean)
#'
#' # check to see whether the count of an internal node is the sum
#' # of counts of its descendant leaves.
#' # here, check the first sample as an example
#'
#' nod <- tx_node(tree = tinyTree, input = rownames(count_tinyTree))
#' d <- cbind.data.frame(node = nod, count = count_tinyTree[, 1])
#'
#' ggtree(tinyTree) %<+% d + geom_text2(aes(label = count))



nodeCount <- function(tipTable, wtree,
                      stree = NULL, fun = sum) {

    if (is.null(tipTable)) {
        stop("tipTable is missing")
    }
    if (!inherits(wtree, "phylo")) {
      stop("wtree: should be a phylo object")
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
        cNode[i, ] <- apply(tipTable[tips.i, ], 2, fun)
    }
    colnames(cNode) <- colnames(tipTable)
    return(rbind(tipTable, cNode))
}
