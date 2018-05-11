#' find descendants (or offsprings)
#'
#' \code{findOS} is to find descendants of an internal node.
#'
#' @param ancestor a
#'
#'
#'
#' @examples
#' data(tinyTree)
#'
#' library(ggtree)
#' ggtree(tinyTree) +
#' geom_text2(aes(label = node)) +
#' geom_hilight(node = 15, fill = 'steelblue', alpha = 0.5)
#'
#' (tips <- findOS(15, tinyTree, only.Tip = TRUE))
#'



findOS <- function(ancestor, tree, only.Tip = TRUE, self.include = TRUE) {
    
    if (length(ancestor) > 1) {
        stop("check whether the argument ancestor is one internal node")
    }
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }
    
    if (!(inherits(ancestor, "character") | inherits(ancestor, "numeric") | inherits(ancestor, 
        "integer"))) {
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
        numA <- tx_node(tree = tree, input = ancestor)
    } else {
        numA <- ancestor
    }
    
    
    # convert to a list.  each element in the list is one path
    desA <- lapply(seq_len(nrow(matN)), FUN = function(x) {
        xx <- match(numA, matN[x, ])
        yy <- ifelse(is.na(xx), 0, xx)
        matN[x, seq_len(yy)]
    })
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
    
    # final output (node number)
    return(res)
    
}

