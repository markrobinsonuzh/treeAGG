#' Transform a phylo object into a matrix.
#'
#' \code{matTree} transforms a phylo tree into a matrix. The entry of the matrix
#' is node number. Each row represents a path connecting a leaf node and the
#' root. The columns are arranged in the order as the path passing the nodes to
#' reach the root.
#'
#' @param tree A phylo object
#' @export
#' @return A matrix
#' @author Ruizhu Huang
#'
#' @examples
#' library(ggtree)
#'
#' ggtree(exTree, branch.length = 'none') +
#'  geom_text2(aes(label = node))
#'
#' data(exTree)
#' # each row of the matrix representing a path.
#' # the first column is leaf nodes; the last non-NA value in a row is the root
#' mat <- matTree(tree = exTree)
#'
matTree <- function(tree) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    # each path connects a tip with the root.
    i <- 1
    mat <- tree$edge
    L1 <- setdiff(mat[, 2], mat[, 1])
    matN <- cbind(L1)
    repeat {
        li <- mat[match(matN[, i], mat[, 2]), 1]
        ll <- length(unique(li))
        if (ll == 1) {
            break
        }
        matN <- cbind(matN, li)
        i <- i + 1
    }
    rownames(matN) <- NULL
    colnames(matN) <- paste("L", seq_len(ncol(matN)), sep = "")

    return(matN)
}

