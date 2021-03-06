#' Calculate true positive rate (TPR) on a tree structure
#'
#' \code{tpr} calculates the true positive rate (TPR) on a tree structure at
#' leaf or node level.
#'
#' @param tree A phylo object
#' @param truth Nodes that have signals (eg. differentally abundant at different
#'   experimental conditions.). If the signals are in different directions (up
#'   and down), then provide the nodes as a list of two members (one up and one
#'   down). One could provide either node number or node label. \strong{Note:}
#'   When TPR at node level is required and only leaf nodes with signal are
#'   provided, then the internal nodes which are shared and only shared by the
#'   provied leaf nodes (signal in the same direction) will be found out and
#'   used with the leaf nodes in the TPR calculation; when TPR at leaf level is
#'   required and the given nodes have internal nodes, then the descendant leaf
#'   nodes will be found out and used in the TPR calculation.
#' @param found Nodes that have been found to have signal (eg. differentally
#'   abundant at different experimental conditions). If the signals are in
#'   different directions (up and down), then provide the nodes as a list of two
#'   members (one up and one down). One could provide either node number or node
#'   label. \strong{Note:} When TPR at node level is required, then the
#'   descendant nodes of the provied nodes (include themselves) will be found
#'   out and used in the TPR calculation; when TPR at leaf level is required,
#'   then the descendant leaf nodes will be found out and used in the TPR
#'   calculation.
#' @param only.Tip A logical value, TRUE or FALSE. If TRUE, true positive rate
#'   is calculated at the leaf (tip) level; otherwise it is calculated at node
#'   level. The default is TRUE
#' @param direction TRUE or FALSE. Default is FALSE. If TRUE, the signal
#'   direction is taken into account; the argument \strong{truth} and
#'   \strong{found} should both be a list of two members and the order of
#'   directions should match.
#'
#' @export
#' @return A true positive rate
#' @author Ruizhu Huang
#' @examples
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
#'             found = c(17, 14), only.Tip = TRUE)
#'  # TPR at the node level
#' tpr2 <- tpr(tree = tinyTree, truth = c(16, 13),
#'             found = c(17, 14), only.Tip = FALSE)
#'
#'
#' # ------------ direction matters-------------------
#' # Found: branches with node 17 (up) and node 14 (down)
#'
#'  # the order of direction for truth and found should be the same
#' (tpr3 <- tpr(tree = tinyTree, truth = list(16, 13),
#'             found = list(16, 13), only.Tip = TRUE,
#'             direction = TRUE))  # correct order
#'
#' (tpr4 <- tpr(tree = tinyTree, truth = list(16, 13),
#'             found = list(13, 16), only.Tip = FALSE,
#'             direction = TRUE)) # wrong order
#'

tpr <- function(tree, truth, found,
                only.Tip = TRUE,
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
                    x <- transNode(tree = tree, input = x,
                                   message = FALSE)
                }
                if (is.character(y)) {
                    y <- transNode(tree = tree, input = y,
                                   message = FALSE)
                }
                .tpr0(tree = tree, truth = x,
                     found = y, only.Tip = only.Tip)
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
            truth <- transNode(tree = tree, input = truth,
                               message = FALSE)
        }
        if (is.character(found)) {
            found <- transNode(tree = tree, input = found,
                               message = FALSE)
        }

        tt <- .tpr0(tree = tree, truth = truth,
                   found = found, only.Tip = only.Tip)
        tpr <- tt[1]/tt[2]
    }

    # return final results
    names(tpr) <- "tpr"
    return(tpr)
}


#' Calculate true positives and positives
#'
#' \code{.tpr0} calculates the number of true positives and the total number of
#' positives on a tree structure at leaf or node level.
#'
#' @param tree A phylo object
#' @param truth Nodes that have signals (eg. differentally abundant at different
#'  experimental conditions.).
#' @param found Nodes that have been found to have signal
#' @param only.Tip A logical value, TRUE or FALSE. If TRUE, true positive rate
#'   is calculated at the leaf (tip) level; otherwise it is calculated at node
#'   level. The default is TRUE
#' @return a vector
#' @keywords internal


.tpr0 <- function(tree,
                 truth = NULL,
                 found = NULL,
                 only.Tip = TRUE) {

    # if no signal (truth = NULL), the true positive is 0 and the positive is 0
    if (is.null(truth)) {
        c(tp = 0, pos = 0)
    } else {
        # if signal exists, check whether the truth has correct input format
        if (!(
            is.character(truth) |
            is.numeric(truth) |
            is.integer(truth)
        )) {
            stop("truth should include character or numeric")
        }

        # if signal exists and has correct input format, but found is null
        # the true positive is 0, and the positive depends on the level
        if (is.null(found)) {
            if (only.Tip) {
                tipT <- findOS(tree = tree, ancestor = truth,
                               only.leaf = TRUE, self.include = TRUE)
                tip <- unlist(tipT)
                c(tp = 0, pos = length(tip))
            } else {
                # the internal node whose descendant leaf nodes all have
                # signal has signal too.
                nodeS <- signalNode(node = truth, tree = tree)
                nodeT <- findOS(tree = tree, ancestor = nodeS,
                                only.leaf = FALSE,
                                self.include = TRUE)
                nod <- unlist(nodeT)
                c(tp = 0, pos = length(nod))
            }
        } else {
            # if signal exists and has correct input format, but found isn't
            # null check whether the input format of found is correct
            if (!(is.character(found) |
                  is.numeric(found) |
                  is.integer(found))) {
                stop("found should include character or numeric")
            }

            # if there is signal and something is found, the true positive and
            # the positive could be calculated after their input formats are
            # correct.
            if (only.Tip) {
                tipT <- findOS(tree = tree, ancestor = truth,
                               only.leaf = TRUE, self.include = TRUE)
                tipF <- findOS(tree = tree, ancestor = found,
                               only.leaf = TRUE, self.include = TRUE)
                # the true positive & positive
                TP <- intersect(unlist(tipT), unlist(tipF))
                tip <- unlist(tipT)
                c(tp = length(TP), pos = length(tip))
            } else {
                # the internal node whose descendant leaf nodes
                # all have signal has signal too.
                nodeS <- signalNode(node = truth, tree = tree)
                # all nodes have signal
                nodeT <- findOS(tree = tree, ancestor = nodeS,
                                only.leaf = FALSE, self.include = TRUE)
                # nodes found with signal
                nodeF <- findOS(tree = tree, ancestor = found,
                                only.leaf = FALSE,
                                self.include = TRUE)
                # true positive & positive
                TP <- intersect(unlist(nodeT), unlist(nodeF))
                nod <- unlist(nodeT)
                c(tp = length(TP), pos = length(nod))
            }
        }
    }
}
