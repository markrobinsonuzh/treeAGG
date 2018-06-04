#' find the node shared by the specified nodes
#'
#' \code{signalNode} is to find nodes which are the ancestors shared by the specified nodes.  The descendant leaves of the found nodes are also the descendant leaves of the specified nodes.
#'
#' @param node a vector of node numbers or node labels
#' @param tree a tree (phylo object)
#' @param label a logical value. If TRUE, the selected node label is output; otherwise the node number is output
#'
#' @return the label of the shared node
#' @export
#'
#' @examples
#'
#' data(tinyTree)
#' library(ggtree)
#'
#' # PLOT tree
#' ggtree(tinyTree,branch.length = 'none')+
#' geom_text2(aes(label = label))
#'
#' ## find the node shared by provided node labels
#' signalNode(node = c('t4','t9'), tree = tinyTree,
#'  label = TRUE)
#'
#' signalNode(node = c('t10','Node_18', 't8'), tree = tinyTree,
#'  label = TRUE)
#'
#' # -----------------------------------------------------
#' ggtree(tinyTree,branch.length = 'none')+
#' geom_text2(aes(label = node))
#'
#' ## find the node shared by provided node numbers
#' signalNode(node = c(2, 3), tree = tinyTree)
#' signalNode(node = c(2, 3, 16), tree = tinyTree)
#'

signalNode <- function(node, tree, label = FALSE) {

  if (!inherits(tree, "phylo")) {
    stop("tree is not a phylo object.")
  }

  if (!is.atomic(node)) {
    stop("node is a vector")
  }

  # transfer node label to node number
  if (inherits(node, "character")) {
    node <- transNode(tree, input = node)
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
  sNode <- rmDesc(node = selN, tree = tree)

  if (label) {
    final <- transNode(tree = tree, input = sNode)
  } else {
    final <- sNode
  }

  return(final)
}
