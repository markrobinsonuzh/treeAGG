#' Convert a tree into multiple paths
#'
#' \code{treePath} convert a phylo object tree into multiple paths. Each starts with a tip, connects all internal nodes, and end with the root.
#'
#' @param tree a phylo object
#'
#' @export
#'
#' @examples
#'
#' library(ape)
#' library(ggtree)
#'
#' # generate a random tree
#' set.seed(1)
#' tree <- rcoal(10)
#' ggtree(tree) + geom_text2(aes(label = node))
#'
#' paths <- pathTree(tree)
#'

pathTree <- function(tree) {
    
    if (is.null(tree) | class(tree) != "phylo") {
        stop("tree: should be a phylo object")
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
    
    # convert to a list.  each element in the list is one path
    listPath <- lapply(seq_len(nrow(matN)), FUN = function(x) {
        y <- matN[x, ]
        y[!is.na(y)]
    })
    return(listPath)
    
}

