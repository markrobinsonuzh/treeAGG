#' prune tree at each node
#'
#' \code{pruneTree} prune the tree at each node.
#'
#' @param wtree an object of phylo class
#'
#' @return a list of phylo object. If the tree ('wtree') provides node and tip labels, the labels are passed to its subtrees; otherwise, the labels are generated automatically with a prefix 'Node_' for internal nodes and 'Tip_' for leaves with the corresponding numbers in the edge component of phylo object. (wtree$edge)
#'
#' @export
#'
#' @examples
#'
#' library(ape)
#' n <- 20
#' tree<- rtree(n)
#' small <- pruneTree(tree, message = TRUE)


pruneTree <- function(tree, message = FALSE) {
    
    if (!inherits(tree, "phylo")) {
        stop("object tree is not of class phylo")
    }
    # paths
    listPath <- pathTree(tree = tree)
    
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
    # --------- prune the tree -----------
    listTree <- lapply(seq_along(nod), FUN = function(x) {
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
        ln.x <- unique(unlist(selPath))
        # select edges
        m.x <- mat[, 1] %in% ln.x
        mat.s <- rank(mat[m.x, ], ties.method = "min")
        mat.x <- matrix(mat.s, nrow = sum(m.x))
        # labels: tip, node
        tip.x <- sort(intersect(mat[m.x, 2], tip))
        nod.x <- sort(intersect(mat[m.x, 1], nod))
        # create phylo object
        list.x <- list(edge = mat.x, tip.label = tipLab[tip.x], edge.length = tree$edge.length[m.x], 
            node.label = nodLab[nod.x - max(tip)], Nnode = length(nod.x))
        class(list.x) <- "phylo"
        list.x
    })
    names(listTree) <- nodLab
    return(listTree)
}


