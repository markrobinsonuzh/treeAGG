#' Transfer between node number and node label
#'
#' \code{transNode} does the transformation between the number and the label of
#' a node on a tree
#'
#' @param tree A phylo object
#' @param input A character or numeric vector representing tree node label(s)
#' or tree node number(s)
#' @param use.original A logical value, TRUE or FALSE. Only required when
#' \strong{input} is numeric. Default is FALSE. There are cases that a tree
#' has no labels for some nodes. If keep.null is TRUE, NA is returned for these
#' nodes; otherwise, the node number (from \strong{input}) will be returned
#' with prefix "Node_" (for internal nodes) or "Leaf_" (for leaf nodes).
#'
#' @export
#' @return a vector
#' @author Ruizhu Huang
#'
#' @examples
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

transNode <- function(tree, input, use.original = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    # node number & tip number
    mat <- tree$edge
    nod <- sort(unique(mat[, 1]))
    tip <- sort(setdiff(mat[, 2], mat[, 1]))

    # check whether the input node number exists in the provided tree
    if (is.numeric(input)) {
        if (!all(input %in% mat)) {
            stop("Node ", input, " can't be found in the ",
                 deparse(substitute(tree)), "\n")
        }
    }

    # tip label
    if (is.null(tree$tip.label)) {
        if(use.original) {
            tipLab <- NULL
        }else{
            tipLab <- paste("leaf_", tip, sep = "")
        }

    } else {
        tipLab <- tree$tip.label
    }
    # node label
    if (is.null(tree$node.label)) {
        if (use.original) {
            nodLab <- NULL
        }else{
            nodLab <- paste("Node_", nod, sep = "")
        }

    } else {
        Labs <- tree$node.label
        if (any(duplicated(Labs))){
            warning("Some internal nodes have same labels")
            nodLab <- paste("Node_", nod, sep = "")
        }else{
            nodLab <- Labs
        }
    }

    comb <- c(tip, nod)
    names(comb) <- c(tipLab, nodLab)

    # transfer from label to number
    if (inherits(input, "character")) {
        if (all(input %in% names(comb))) {
            final <- comb[input]
        } else {
            stop("The nodes ", paste(input[!input %in% names(comb)], collapse = ", "),
                 " could not be found in the tree. \n Node numbers or Node labels are
           required but not a mixture of both")
        }

    } else {
        final <- names(comb[match(input, comb)])
    }

    return(final)

}
