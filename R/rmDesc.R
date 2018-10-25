#' Remove nodes which are the descendants of other nodes
#'
#' \code{rmDesc} removes branches, which are the subbranches of the others.
#'
#' @param tree A phylo object.
#' @param node A vector of node labels or node numbers
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"Node_"} to the node number if the node is
#'   an internal node or adding a prefix \code{"Leaf_"} if the node is a leaf
#'   node.
#'
#' @export
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#' @author Ruizhu Huang
#' @examples {
#' data(tinyTree)
#' library(ggtree)
#' # PLOT
#' ggtree(tinyTree, branch.length = 'none') +
#'  geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'  geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#'
#' # find the shared nodes from the tree plot
#' aV <- c(5, 4, 18, 3, 2)
#' # final result
#' (rn <- rmDesc(tree = tinyTree, node = aV))
#'}
#'

rmDesc <- function(tree, node, use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (is.character(node)) {
        node <- transNode(tree = tree, input = node,
                          message = FALSE)
    } else {
        node <- node
    }

    # find descendant leaves
    ListTip <- lapply(as.list(node),
                      FUN = function(x) {
                          findOS(ancestor = x, tree = tree,
                                 only.Tip = TRUE, self.include = TRUE)
    })

    # remove duplicates
    ListTip1 <- ListTip[!duplicated(ListTip)]

    # find out elements that are not the subset of the others
    ind <- rep(TRUE, length(ListTip1))
    for (i in seq_along(ListTip1)) {
        xi <- ListTip1[[i]]
        ti <- lapply(ListTip1, FUN = function(x) {
            all(xi %in% x)
        })
        ti <- unlist(ti)
        if (sum(ti) > 1) {
            ind[i] <- FALSE
        }
    }
    out <- node[!duplicated(ListTip)][ind]

    # return a vector of the found node (the node number of the node)
    # name the vector with the node label
    names(out) <- transNode(tree = tree, input = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)

}
