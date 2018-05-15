#' find the ancestors of specified nodes
#'
#' \code{findAncestor} is to find the ancestor in the n generation before specified nodes.
#'
#' @param tree a phylo object
#' @param node a vector of node numbers or node labels
#' @param level a vector of numbers to define n generation before the specified nodes
#' @param treeMat a matrix with each row representing a path. The entry is node number. The first column is leave node numbers. The columns are organized as the order of nodes in the paths connecting leaves and the root. Default is null, treeMat is generated automatically. The treeMat is provided when the tree is very large to save running time.
#'
#'@export
#'@return a vector of node numbers
#'
#'@examples
#' library(ggtree)
#' data(exTree)
#' ggtree(exTree, branch.length = 'none') +
#'  geom_text2(aes(label = node))
#'
#'  findAncestor(tree = exTree, node = c(53, 61), level = 1)


findAncestor <- function(tree, node, level, treeMat = NULL) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (is.null(treeMat)) {
        treeMat <- matTree(tree)
    } else {
        treeMat <- treeMat
    }

    if (inherits(node, "character")) {
        aggNod <- tx_node(tree = tree, input = node)
    } else {
        aggNod <- node
    }

    if (length(level) == 1) {
        level <- rep(level, length(node))
    } else {
        if (length(level) == length(node)) {
            level <- level
        } else {
            stop("the length of level is not equal to the length of node")
        }
    }

    selNod <- lapply(seq_along(aggNod), FUN = function(x) {
        # find rows with selected nodes
        ind <- which(treeMat == aggNod[x], arr.ind = TRUE)
        valP <- apply(treeMat, 1, FUN = function(x) {
            max(which(!is.na(x)))
        })
        selP <- valP[ind[, 1]]

        # move up levels as specified until the root
        vv <- ind[, 2] + level[x]
        if (all(vv < selP)) {
            ind.f <- cbind(ind[, 1], vv)
        } else {
            # v2 <- ifelse(vv <= selP, vv, selP)
            stop("The level specified for nodes ", node[x], " exceed the root level; try a small level value. \n")
            # ind.f <- cbind(ind[, 1], v2)
        }

        # there is only one path to go to the root when the starting point is fixed
        vu <- unique(treeMat[ind.f])
        if (length(vu) > 1) {
            stop("more than one node has been found.")
        } else {
            return(vu)
        }
    })

    final <- unlist(selNod)
    return(final)
}
