#' Find the ancestors of specified nodes
#'
#' \code{findAncestor} finds the ancestor in the nth generation above
#' specified nodes.
#'
#' @param tree A phylo object
#' @param node A vector of node numbers or node labels
#' @param level A vector of numbers to define nth generation before the
#' specified nodes
#' @param return "number" (return the node number) or "label" (return the node
#'   label).
#' @param use.alias A logical value, TRUE or FALSE. This is an optional argument
#'   that only requried when \code{return = "label"}. The default is FALSE, and
#'   the node label would be returned; otherwise, the alias of node label would
#'   be output. The alias of node label is created by adding a prefix
#'   \code{"Node_"} to the node number if the node is an internal node or adding
#'   a prefix \code{"Leaf_"} if the node is a leaf node.
#' @export
#' @return a vector of node numbers
#' @author Ruizhu Huang
#'
#' @examples
#' library(ggtree)
#' data(exTree)
#' ggtree(exTree, branch.length = 'none') %>%
#'     scaleClade(node = 52, scale = 10)+
#'  geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'  geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#'  findAncestor(tree = exTree, node = c(53, 61), level = 1)

findAncestor <- function(tree, node, level,
                         return = c("number", "label"),
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

    # final output (node number or label)
    return <- match.arg(return)
    switch(return,
           number = final,
           label = transNode(tree = tree, input = final,
                             use.alias = use.alias,
                             message = FALSE))

}
