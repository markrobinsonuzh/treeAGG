#' Calculate counts at internal nodes
#'
#' \code{nodeCount} calculates the counts at internal nodes. The count of an
#' internal node is the sum of counts at its descendant leaves.
#'
#' @param tree A phylo object
#' @param data A matrix or data frame. A count table from real data.
#' @param fun A function to create the count of an internal node based on the
#'   counts at its descendant leaf nodes. The default is sum
#'
#' @importFrom utils head
#' @export
#' @return A count table (matrix class) with a row representing a node and a
#'   column representing a sample.
#' @author Ruizhu Huang
#'
#'
#' @examples
#' data("tinyTree")
#' library(ggtree)
#' (p <- ggtree(tinyTree) + geom_text(aes(label = label)))
#'
#' count <- matrix(rpois(100, 10), nrow =10)
#' rownames(count) <- tinyTree$tip.label
#' colnames(count) <- paste(c("C1_", "C2_"),
#' c(1:5, 1:5), sep = "")
#'
#' count_tinyTree <- nodeCount(data = count,
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

nodeCount <- function(tree, data, fun = sum) {
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
