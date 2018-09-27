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
#' @param return "number" (return the node number) or "label" (return the node
#'   label).
#'
#' @export
#' @return A vector
#' @author Ruizhu Huang
#'
#' @examples
#' data(tinyTree)
#'
#' library(ggtree)
#' ggtree(tinyTree) +
#' geom_text2(aes(label = node)) +
#' geom_hilight(node = 15, fill = 'steelblue', alpha = 0.5)
#'
#' (tips <- findOS(tree = tinyTree, ancestor = 15, only.Tip = TRUE))

findOS <- function(tree,
                   ancestor,
                   only.Tip = TRUE,
                   self.include = TRUE,
                   return = c("number", "label")) {
    if (length(ancestor) > 1) {
        stop("ancestor should have length equal to one")
    }
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (!(
        inherits(ancestor, "character") |
        inherits(ancestor, "numeric") | inherits(ancestor,
                                                 "integer")
    )) {
        stop("ancestor should be character or numeric")
    }
    # each path connects a tip with the root.
    i <- 1
    mat <- tree$edge
    L1 <- setdiff(mat[, 2], mat[, 1])
    matN <- cbind(L1)
    repeat {
        li <- mat[match(matN[, i], mat[, 2]), 1]
        ll <- length(unique(li))
        if (ll == 1) {
            break
        }
        matN <- cbind(matN, li)
        i <- i + 1
    }
    rownames(matN) <- colnames(matN) <- NULL

    if (inherits(ancestor, "character")) {
        numA <- transNode(tree = tree, input = ancestor,
                          use.original = FALSE,
                          message = FALSE)
    } else {
        numA <- ancestor
        if (!numA %in% mat) {
            stop("Node ",
                 numA,
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
    return <- match.arg(return)
    switch(return,
           number = res,
           label = transNode(tree = tree, input = res,
                             use.original = FALSE,
                             message = FALSE))

}
