#' Calculate false discovery rate (fdr) on a tree structure
#'
#' \code{fdr} calculates the false discovery rate (fdr) on a tree structure at
#' the leaf or node level.
#'
#' @param tree A phylo object
#' @param truth Nodes that have signals (eg. differentally abundant at different
#'   experimental conditions.). If the signals are in different directions (up
#'   and down), then provide the nodes as a list of two members (one up and one
#'   down). One could provide either node number or node label. \strong{Note:}
#'   When fdr at node level is required and only leaf nodes with signal are
#'   provided, then the internal nodes which are shared and only shared by the
#'   provied leaf nodes (signal in the same direction) will be found out and
#'   used with the leaf nodes in the fdr calculation; when fdr at leaf level is
#'   required and the given nodes have internal nodes, then the descendant leaf
#'   nodes will be found out and used in the fdr calculation.
#' @param found Nodes that have been found to have signal (eg. differentally
#'   abundant at different experimental conditions). If the signals are in
#'   different directions (up and down), then provide the nodes as a list of two
#'   members (one up and one down). One could provide either node number or node
#'   label. \strong{Note:} When fdr at node level is required, then the
#'   descendant nodes of the provied nodes (include themselves) will be found
#'   out and used in the fdr calculation; when fdr at leaf level is required,
#'   then the descendant leaf nodes will be found out and used in the fdr
#'   calculation.
#' @param level If "leaf", false discovery rate is calculated at leaf level; if
#'   "node", it is calculated at node level.
#' @param direction TRUE or FALSE. Default is FALSE. If TRUE, the signal
#'   direction is taken into account; the argument \strong{truth} and
#'   \strong{found} should both be a list of two members and the order of
#'   directions should match.
#'
#' @export
#' @return A false discovery rate
#'
#' @examples
#'
#' library(ggtree)
#' data("tinyTree")
#' ggtree(tinyTree) + geom_text2(aes(label = node))
#'
#' # Truth: two branches have differential
#' # abundance under different conditions.
#' # branch with node 16 (direction: increase)
#' # branch with node 13 (direction: decrease)
#'
#'
#' # ------------ ignore direction -------------------
#'  # Found: branches with node 17 and node 14
#'  # fdr at the tip level
#' fdr1 <- fdr(tree = tinyTree, truth = c(16, 13),
#'             found = c(17, 14), level = "leaf")
#'  # fdr at the node level
#' fdr2 <- fdr(tree = tinyTree, truth = c(16, 13),
#'             found = c(17, 14), level = "node")
#'
#'
#' # ------------ direction matters-------------------
#' # Found: branches with node 17 (up) and node 14 (down)
#'
#'  # the order of direction for truth and found should be the same
#' (fdr3 <- fdr(tree = tinyTree, truth = list(16, 13),
#'             found = list(16, 13), level = "leaf",
#'             direction = TRUE))  # correct order
#'
#' (fdr4 <- fdr(tree = tinyTree, truth = list(16, 13),
#'             found = list(13, 16), level = "node",
#'             direction = TRUE)) # wrong order
#'
#'
#'

fdr <- function(tree, truth, found,
                level = c("node", "leaf"),
                direction = FALSE) {

  if (!inherits(tree, "phylo")) {
    stop("tree: should be a phylo object")
  }

  # if signal direction is taken into account.
  if (direction) {
    # it requires list input for both truth and found
    # the length of list should equal to 2 (direction up & down)
    if (inherits(truth, "list") &&
        inherits(found, "list") &&
        length(truth) == length(found)) {

      tt <- mapply(function(x, y) {
        # transfer node label to node number
        if (is.character(x)) {
          x <- transNode(tree = tree, input = x)
        }
        if (is.character(y)) {
          y <- transNode(tree = tree, input = y)
        }

        fdr0(tree = tree, truth = x,
             found = y, level = level)
      }, x = truth, y = found)

      fdr <- rowSums(tt)[1]/rowSums(tt)[2]
    } else {
      stop("check: \n
           truth is a list of two members? \n
           found is a list of two members? \n ")
    }

    # if signal direction isn't taken into account.
    } else {
      # transfer node label to node number
      if (is.character(truth)) {
        truth <- transNode(tree = tree, input = truth)
      }
      if (is.character(found)) {
        found <- transNode(tree = tree, input = found)
      }
      tt <- fdr0(tree = tree, truth = truth,
                 found = found, level = level)
      fdr <- tt[1]/tt[2]
    }

  # return final results
  names(fdr) <- "fdr"
  return(fdr)
}


#' Calculate the number of false discoveries and the total number of discoveries on a tree structure
#'
#' \code{fdr0} calculates the number of false discoveries and the total number of discoveries at leaf or node level on a tree structure .
#'
#' @param tree A phylo object
#' @param truth Nodes that have signals (eg. differentally abundant at different experimental conditions.).
#' @param found Nodes that have been found to have signal
#' @param level If "leaf", false discovery rate is calculated at leaf level; if "node", it is calculated at node level.
#' @return

fdr0 <- function(tree,
                 truth = NULL,
                 found = NULL,
                 level = c("node", "leaf")) {

  # if no discovery (found = NULL), the false discovery is 0 and the discovery is 0
  if (is.null(found)) {
    c(fd = 0, disc = 0)
  } else {
    # if discovery exists, check whether the discovery has correct input format
    if (!(
      inherits(found, "character") |
      inherits(found, "numeric") |
      inherits(found, "integer")
    )) {
      stop("found should include character or numeric")
    }

    # if discovery exists and has correct input format, but truth is null
    # then false discovery equals to discovery and the number depends on the level
    if (is.null(truth)) {
      level <- match.arg(level)
      switch(level,
             leaf = {
               tipF <- lapply(
                 found,
                 findOS,
                 tree = tree,
                 only.Tip = TRUE,
                 self.include = TRUE
               )
               tip <- unlist(tipF)
               c(fp = length(tip), disc = length(tip))
             },
             node = {
               # the internal node whose descendant leaf nodes all have signal
               # has signal too.
               nodeF <- lapply(
                 found,
                 findOS,
                 tree = tree,
                 only.Tip = FALSE,
                 self.include = TRUE
               )
               nod <- unlist(nodeF)
               c(fp = length(nod), disc = length(nod))
             })
    } else {

      # if discovery (found) exists and has correct input format, and truth isn't null
      # check whether the input format of truth is correct
      if (!(
        inherits(truth, "character") |
        inherits(truth, "numeric") |
        inherits(truth, "integer")
      )) {
        stop("truth should include character or numeric")
      }

      # if both found and truth are not NULL, the false discovery and the        discovery could be calculated after their input formats are correct as below.
      level <- match.arg(level)
      switch(level,
             leaf = {
               tipT <- lapply(
                 truth,
                 findOS,
                 tree = tree,
                 only.Tip = TRUE,
                 self.include = TRUE
               )
               tipF <- lapply(
                 found,
                 findOS,
                 tree = tree,
                 only.Tip = TRUE,
                 self.include = TRUE
               )
               # the true positive & positive
               fp <- setdiff(unlist(tipF), unlist(tipT))
               tip <- unlist(tipF)
               c(fp = length(fp), disc = length(tip))
             },
             node = {
               # the internal node whose descendant leaf nodes
               # all have signal has signal too.
               nodeS <- signalNode(node = truth,
                                   tree = tree,
                                   label = TRUE)
               # all nodes have signal
               nodeT <- lapply(
                 nodeS,
                 findOS,
                 tree = tree,
                 only.Tip = FALSE,
                 self.include = TRUE
               )
               # nodes found with signal
               nodeF <- lapply(
                 found,
                 findOS,
                 tree = tree,
                 only.Tip = FALSE,
                 self.include = TRUE
               )
               # false discovery & discovery
               fp <- setdiff(unlist(nodeF), unlist(nodeT))
               nod <- unlist(nodeF)
               c(fp = length(fp), disc = length(nod))
             })
    }
  }
}
