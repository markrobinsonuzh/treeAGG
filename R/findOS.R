#' Find descendants (or offsprings)
#'
#' \code{findOS} finds descendants of a node.
#'
#' @param ancestor An internal node. It could be the node number or the node
#'   label.
#' @param tree A phylo object.
#' @param only.Tip A logical value, TRUE or FALSE. The default is TRUE. If
#'   default, only the leaf nodes in the descendant nodes would be returned.
#' @param self.include A logical value, TRUE or FALSE. The default is TRUE. If
#'   default, the node specified in \strong{ancestor} is included. The leaf node
#'   itself is returned as its descendant.
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"Node_"} to the node number if the node is
#'   an internal node or adding a prefix \code{"Leaf_"} if the node is a leaf
#'   node.
#' @export
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#' @author Ruizhu Huang
#'
#' @examples
#' data(tinyTree)
#'
#' library(ggtree)
#' ggtree(tinyTree) +
#' geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7) +
#' geom_hilight(node = 17, fill = 'steelblue', alpha = 0.5) +
#' geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7)
#'
#' (tips <- findOS(tree = tinyTree, ancestor = 17, only.Tip = TRUE))

findOS <- function(tree,
                   ancestor,
                   only.Tip = TRUE,
                   self.include = TRUE,
                   use.alias = FALSE) {
    if (length(ancestor) > 1) {
        stop("ancestor should have length equal to one")
    }
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (!(is.character(ancestor) |
          is.numeric(ancestor) |
          is.integer(ancestor))) {
        stop("ancestor should be character or numeric")
    }
    # the edge matrix
    mat <- tree$edge
    matN <- matTree(tree = tree)

    if (is.character(ancestor)) {
        numA <- transNode(tree = tree, input = ancestor,
                          use.alias = TRUE,
                          message = FALSE)
    } else {
        numA <- ancestor
        if (!numA %in% mat) {
            stop("Node ", numA,
                 " can't be found in the ",
                 deparse(substitute(tree)),
                 "\n")
        }

    }

    # convert to a list. (each element in the list is one path)
    desA <- lapply(
        seq_len(nrow(matN)),
        FUN = function(x) {
            xx <- match(numA, matN[x, ])
            yy <- ifelse(is.na(xx), 0, xx)
            matN[x, seq_len(yy)]
        }
    )
    # descendants: internal nodes & tips
    desA <- unique(unlist(desA))

    # descendants: tips
    tipA <- unique(setdiff(desA, mat[, 1]))

    res <- if (self.include) {
        if (only.Tip) {
            tipA
        } else {
            desA
        }
    } else {
        if (only.Tip) {
            setdiff(tipA, numA)
        } else {
            setdiff(desA, numA)
        }
    }

    # final output (node number or label)
    out <- res
    names(out) <- transNode(tree = tree, input = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
}
