#' find the optimal nodes to short result.
#'
#' \code{signalNode} is to represent some nodes with their ancestor to make
#' result as short as possible. The ancestors share exactly the same leaves as
#' the original nodes.
#'
#' @param tree A tree (phylo object)
#' @param node A vector of node numbers or node labels
#' @param return "label" or "number". Default is "number", the node number is
#'   returned. If "label", the selected node label is returned.
#' @param use.alias A logical value, TRUE or FALSE. This is an optional argument
#'   that only requried when \code{return = "label"}. The default is FALSE, and
#'   the node label would be returned; otherwise, the alias of node label would
#'   be output. The alias of node label is created by adding a prefix
#'   \code{"Node_"} to the node number if the node is an internal node or adding
#'   a prefix \code{"Leaf_"} if the node is a leaf node.
#' @export
#' @return The label of the shared node
#' @author Ruizhu Huang
#' @examples
#'
#' data(tinyTree)
#' library(ggtree)
#'
#' # PLOT tree
#' # The node labels are in orange texts and the node numbers are in blue
#' ggtree(tinyTree,branch.length = 'none')+
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#' ## find the node shared by provided node labels
#' signalNode(node = c('t4','t9'), tree = tinyTree,
#'  return = "label")
#' signalNode(node = c('t4','t9'), tree = tinyTree,
#'  return = "number")
#' signalNode(node = c('t10','Node_18', 't8'), tree = tinyTree,
#'  return = "label", use.alias = FALSE)
#' signalNode(node = c('t10','Node_18', 't8'), tree = tinyTree,
#'  return = "label", use.alias = TRUE)
#'
#' ## find the node shared by provided node numbers
#' signalNode(node = c(2, 3), tree = tinyTree)
#' signalNode(node = c(2, 3, 16), tree = tinyTree)
#'

signalNode <- function(tree, node, return = c("number", "label"),
                       use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree is not a phylo object.")
    }

    if (!is.atomic(node)) {
        stop("node is a vector")
    }

    # transfer node label to node number
    if (inherits(node, "character")) {
        node <- transNode(tree, input = node,
                          message = FALSE)
    } else {
        node <- node
    }

    # path matrix
    mat <- matTree(tree)

    if (!all(node %in% as.vector(mat))) {
        stop("Some nodes could not be found in the tree")
    }

    # select paths which include the input nodes
    ind <- apply(mat, 1, FUN = function(x) {
        any(x %in% node)
    })
    # select nodes which only exist in the selected paths
    selN <- setdiff(as.vector(mat[ind, ]), as.vector(mat[!ind, ]))
    # remove nodes which are descendants of any others
    matI <- mat[ind,, drop = FALSE]
    selF <- apply(matI, MARGIN = 1,
                FUN = function(x){
                    y <- x %in% selN
                    x[sum(y)]
                })
    sNode <- unique(selF)
    sNode <- sNode[!is.na(sNode)]

    # final output (node number or label)
    return <- match.arg(return)
    switch(return,
           number = sNode,
           label = transNode(tree = tree, input = sNode,
                             use.alias = use.alias,
                             message = FALSE))
}
