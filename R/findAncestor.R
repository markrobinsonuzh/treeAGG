#' Find the ancestors of specified nodes
#'
#' \code{findAncestor} finds the ancestor in the nth generation above
#' specified nodes.
#'
#' @param tree A phylo object
#' @param node A vector of node numbers or node labels
#' @param level A vector of numbers to define nth generation before the
#' specified nodes
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"Node_"} to the node number if the node is
#'   an internal node or adding a prefix \code{"Leaf_"} if the node is a leaf
#'   node.
#' @export
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#' @author Ruizhu Huang
#'
#' @examples
#' library(ggtree)
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = 'none') +
#'  geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'  geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#'  findAncestor(tree = tinyTree, node = c(18, 13), level = 1)

findAncestor <- function(tree, node, level,
                         use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    # convert a tree to a matrix
    # each row is a path connecting the root and a leaf
    treeMat <- matTree(tree)


    if (is.character(node)) {
        aggNod <- transNode(tree = tree, input = node,
                            use.alias = TRUE,
                            message = FALSE)
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
        selP <- valP[ind[, "row"]]

        # move up levels as specified until the root
        vv <- ind[, "col"] + level[x]
        if (all(vv <= selP)) {
            ind.f <- cbind(ind[, "row"], vv)
        } else {
            stop("The level specified for nodes ", node[x], " exceed the root
                 level; try a small level value. \n")

        }

        # there is only one path to go to the root when the starting point is
        # fixed
        vu <- unique(treeMat[ind.f])
        if (length(vu) > 1) {
            stop("more than one node has been found.")
        } else {
            return(vu)
        }
    })

    final <- unlist(selNod)

    # return a vector of the found node (the node number of the node)
    # name the vector with the node label
    out <- final
    names(out) <- transNode(tree = tree, input = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
}
