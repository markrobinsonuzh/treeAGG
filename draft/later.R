#' prune tree at each node
#'
#' \code{pruneTree} prune the tree at each node.
#'
#' @param wtree an object of phylo class
#'
#' @return a list of phylo object
#'
#' @export
#'
#' @examples
#'
#' library(ape)
#' n <- 20
#' Lab <- paste('Node',1:(n-1),sep='')
#'
#' # entire tree
#' Tree <- TreeW<- rtree(n)
#'
#' # Tree1 has node labels
#' TreeW$node.label <- Lab
#'
#'
#' small.w <- pruneTree(Tree1, label = TRUE)
#' small <- pruneTree(Tree, label = FALSE)
#' small2 <- pruneTree(Tree, label = TRUE)

pruneTree <- function(tree, message = FALSE, label = FALSE) {

  if (!inherits(tree, "phylo")) {
    stop("object tree is not of class phylo")
  }
  # paths
  listPath <- pathTree(tree = tree)

  # node number & tip number
  mat <- tree$edge
  nod <- unique(mat[, 1])
  tip <- setdiff(mat[, 2], mat[, 1])

  # --------- match node number and their labels ----------- node label
  nodLab <- paste("Node_", nod, sep = "")
  if (is.null(tree$node.label)) {
    names(nod) <- nodLab
  } else {
    names(nod) <- tree$node.label
  }
  # tip label
  tipLab <- paste("Tip_", tip, sep = "")
  if (is.null(tree$tip.label)) {
    names(tip) <- tipLab
  } else {
    names(tip) <- tree$tip.label
  }
  # combine label
  Lab <- c(nod, tip)

  # -----------decide to return nodes or labels ------- list of subtrees (output node
  # number)
  listNode <- lapply(seq_along(nod), FUN = function(x) {
    if (message) {
      cat(x, "out of", length(nod), "nodes finished \n")
    }
    nod.x <- nod[x]
    selPath <- lapply(listPath, FUN = function(y) {
      if (nod.x %in% y) {
        y[seq_len(which(y == nod.x))]
      } else {
        NULL
      }
    })
    unique(unlist(selPath))
  })

  if (label) {
    # list of subtrees (output node label)
    listLab <- lapply(seq_along(listNode), FUN = function(x) {
      list.x <- listNode[[x]]
      lab.x <- Lab[match(list.x, Lab)]
      names(lab.x)
    })

    return(listLab)
  } else {
    return(listNode)
  }
}






# paths
listPath <- pathTree(tree = tree)

# node number & tip number
mat <- tree$edge
nod <- sort(unique(mat[, 1]))
tip <- sort(setdiff(mat[, 2], mat[, 1]))

# --------- match node number and their labels -----------
listNode <- lapply(seq_along(nod), FUN = function(x) {
  if (message) {
    cat(x, "out of", length(nod), "nodes finished \n")
  }
  nod.x <- nod[x]
  selPath <- lapply(listPath, FUN = function(y) {
    if (nod.x %in% y) {
      y[seq_len(which(y == nod.x))]
    } else {
      NULL
    }
  })
  unique(unlist(selPath))
})
