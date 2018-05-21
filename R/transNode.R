#' transfer between node number and node label
#'
#' \code{transNode} is to do the transfermation between the number and the label of a node on a tree
#'
#' @param tree a phylo object
#' @param input a tree node label or a tree node number
#'
#' @export
#'
#' @examples
#'
#'
#' library(ggtree)
#'
#' data(tinyTree)
#'
#' ggtree(tinyTree, branch.length = 'none') +
#' geom_text2(aes(label = label), hjust = -0.3) +
#' geom_text2(aes(label = node), vjust = -0.8,
#' hjust = -0.3, color = 'blue')
#'
#' #check whether the node number and node label are matched
#' transNode(tinyTree, input = c(11, 2, 4))
#'
#' transNode(tinyTree, input = c('Node_16', 'Node_11'))
#'
#'

transNode <- function(tree, input) {

  if (!inherits(tree, "phylo")) {
    stop("tree: should be a phylo object")
  }

  # node number & tip number
  mat <- tree$edge
  nod <- sort(unique(mat[, 1]))
  tip <- sort(setdiff(mat[, 2], mat[, 1]))

  # tip label
  if (is.null(tree$tip.label)) {
    tipLab <- paste("Tip_", tip, sep = "")
  } else {
    tipLab <- tree$tip.label
  }
  # node label
  if (is.null(tree$node.label)) {
    nodLab <- paste("Node_", nod, sep = "")
  } else {
    nodLab <- tree$node.label
  }

  comb <- c(tip, nod)
  names(comb) <- c(tipLab, nodLab)

  # transfer from label to number
  if (inherits(input, "character")) {
    if (all(input %in% names(comb))) {
      final <- comb[input]
    } else {
      stop("The nodes ", paste(input[!input %in% names(comb)], collapse = ", "),
           " could not be found in the tree. \n Node numbers or Node labels are required but not a mixture of both")
    }

  } else {
    final <- names(comb[match(input, comb)])
  }

  return(final)

}
