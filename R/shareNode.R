#' Find the share node
#'
#' \code{shareNode} is to find the node where the specified nodes
#' first meet.
#'
#' @param tree A  phylo object.
#' @param node A vector of node numbers or node labels.
#' @param label A logical value. If TRUE, the selected node label is output;
#'   otherwise the node number is output
#'
#' @return The label of the shared node
#' @export
#'
#' @examples
#'
#' data(tinyTree)
#' library(ggtree)
#' # PLOT tree
#' ggtree(tinyTree,branch.length = 'none')+
#' geom_text2(aes(label = label))
#'
#' ## find the node shared by provided node labels
#' shareNode(node = c('t4','t9'), tree = tinyTree,
#'  label = TRUE)
#'
#' shareNode(node = c('t10','Node_17'), tree = tinyTree,
#'  label = TRUE)
#'
#' # -----------------------------------------------------
#' ggtree(tinyTree,branch.length = 'none')+
#' geom_text2(aes(label = node))
#'
#' ## find the node shared by provided node numbers
#' shareNode(node = c(2, 3), tree = tinyTree)

shareNode <- function(tree, node, label = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree is not a phylo object.")
    }

    if (!is.atomic(node)) {
        stop("node is a vector")
    }

    # transfer node label to node number
    if (inherits(node, "character")) {
        node <- transNode(tree, input = node)
    } else {
        node <- node
    }

    # path matrix
    mat <- matTree(tree)
    ind <- apply(mat, 1, FUN = function(x) {
        any(x %in% node)
    })
    matN <- mat[ind, ]
    path <- lapply(seq_len(nrow(matN)), FUN = function(x) {
        xx <- matN[x, ]
        xx[!is.na(xx)]
    })
    # ancestors
    loc <- Reduce(intersect, path)

    # the ancestor on the lowest level (the root has the highest level)
    vec <- as.vector(mat)
    count <- table(vec, useNA = "no")[loc]
    df <- as.data.frame(count, stringsAsFactors = FALSE)

    # select the node with the lowest frequency.  closest to the leaf level.
    sNode <- as.numeric(df[df$Freq == min(df$Freq), 1])

    if (label) {
        final <- transNode(tree = tree, input = sNode)
    } else {
        final <- sNode
    }

    return(final)
}
