#' Calculate true positive rate (TPR) on a tree structure
#'
#' \code{tpr} calculates the true positive rate (TPR) on a tree structure at leaf or node level.
#'
#' @param tree A phylo object
#' @param truth Nodes that have signals (eg. differentally abundant at different experimental conditions.). If the signals are in different directions (up and down), then provide the nodes as a list of two members (one up and one down). One could provide either node number or node label. \strong{Note:} When TPR at node level is required and only leaf nodes with signal are provided, then the internal nodes which are shared and only shared by the provied leaf nodes (signal in the same direction) will be found out and used with the leaf nodes in the TPR calculation; when TPR at leaf level is required and the given nodes have internal nodes, then the descendant leaf nodes will be found out and used in the TPR calculation.
#' @param found Nodes that have been found to have signal (eg. differentally abundant at different experimental conditions). If the signals are in different directions (up and down), then provide the nodes as a list of two members (one up and one down). One could provide either node number or node label. \strong{Note:} When TPR at node level is required, then the descendant nodes of the provied nodes (include themselves) will be found out and used in the TPR calculation; when TPR at leaf level is required, then the descendant leaf nodes will be found out and used in the TPR calculation.
#' @param level If "leaf", true positive rate is calculated at leaf level; if "node", it is calculated at node level.
#' @param direction TRUE or FALSE. Default is FALSE. If TRUE, the signal direction is taken into account; the argument \strong{truth} and \strong{found} should both be a list of two members and the order of directions should match.
#'
#' @export
#' @return A true positive rate
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
#' # ------------ ignore direction -------------------
#' # Found: branches with node 17 and node 14
#'  # TPR at the tip level
#' tpr1 <- tpr(tree = tinyTree, truth = c(16, 13),
#'             found = c(17, 14), level = "leaf")
#'  # TPR at the node level
#' tpr2 <- tpr(tree = tinyTree, truth = c(16, 13),
#'             found = c(17, 14), level = "node")
#'
#'
#' # ------------ direction matters-------------------
#' # Found: branches with node 17 (up) and node 14 (down)
#'
#'  # the order of direction for truth and found should be the same
#' (tpr3 <- tpr(tree = tinyTree, truth = list(16, 13),
#'             found = list(16, 13), level = "leaf",
#'             direction = TRUE))  # correct order
#'
#' (tpr4 <- tpr(tree = tinyTree, truth = list(16, 13),
#'             found = list(13, 16), level = "node",
#'             direction = TRUE)) # wrong order
#'

tpr <- function(tree, truth, found,
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
        if (is.character(x)) {
          x <- transNode(tree = tree, input = x)
        }
        if (is.character(y)) {
          y <- transNode(tree = tree, input = y)
        }
        tpr0(tree = tree, truth = x,
             found = y, level = level)
      }, x = truth, y = found)
      
      tpr <- rowSums(tt)[1]/rowSums(tt)[2]
    } else {
      stop("check: \n
            truth is a list of two members? \n
            found is a list of two members? \n ")
    }
    
  } else {
    # if signal direction isn't taken into account.
    
    if (is.character(truth)) {
      truth <- transNode(tree = tree, input = truth)
    }
    if (is.character(found)) {
      found <- transNode(tree = tree, input = found)
    }

    tt <- tpr0(tree = tree, truth = truth,
               found = found, level = level)
    tpr <- tt[1]/tt[2]
  }

  # return final results
  names(tpr) <- "tpr"
  return(tpr)
}


#' Calculate the number of true positives and the total number of positives on a tree structure
#'
#' \code{tpr0} calculates the number of true positives and the total number of positives on a tree structure at leaf or node level.
#'
#' @param tree A phylo object
#' @param truth Nodes that have signals (eg. differentally abundant at different experimental conditions.).
#' @param found Nodes that have been found to have signal
#' @param level If "leaf", true positive rate is calculated at leaf level; if "node", it is calculated at node level.
#'

tpr0 <- function(tree,
                 truth = NULL,
                 found = NULL,
                 level = c("node", "leaf")) {

  # if no signal (truth = NULL), the true positive is 0 and the positive is 0
  if (is.null(truth)) {
    c(tp = 0, pos = 0)
  } else {
    # if signal exists, check whether the truth has correct input format
    if (!(
      inherits(truth, "character") |
      inherits(truth, "numeric") |
      inherits(truth, "integer")
    )) {
      stop("truth should include character or numeric")
    }
    
    # if signal exists and has correct input format, but found is null
    # the true positive is 0, and the positive depends on the level
    if (is.null(found)) {
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
               tip <- unlist(tipT)
               c(tp = 0, pos = length(tip))
             },
             node = {
               # the internal node whose descendant leaf nodes all have signal
               # has signal too.
               nodeS <- signalNode(node = truth,
                                   tree = tree,
                                   label = TRUE)
               nodeT <- lapply(
                 nodeS,
                 findOS,
                 tree = tree,
                 only.Tip = FALSE,
                 self.include = TRUE
               )
               nod <- unlist(nodeT)
               c(tp = 0, pos = length(nod))
             })
    } else {
      # if signal exists and has correct input format, but found isn't null
      # check whether the input format of found is correct
      if (!(
        inherits(found, "character") |
        inherits(found, "numeric") | inherits(found,
                                              "integer")
      )) {
        stop("found should include character or numeric")
      }

      # if there is signal and also something found, the true positive and the positive could be calculated after their input formats are correct.
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
               TP <- intersect(unlist(tipT), unlist(tipF))
               tip <- unlist(tipT)
               c(tp = length(TP), pos = length(tip))
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
               # true positive & positive
               TP <- intersect(unlist(nodeT), unlist(nodeF))
               nod <- unlist(nodeT)
               c(tp = length(TP), pos = length(nod))
             })
    }
  }
}
