#' Find the share node
#'
#' \code{shareNode} is to find the node where the specified nodes first meet.
#'
#' @param tree A  phylo object.
#' @param node A vector of node numbers or node labels.
#' @param return "number" (return the node number) or "label" (return the node
#'   label).
#' @param use.alias A logical value, TRUE or FALSE. This is an optional argument
#'   that only requried when \code{return = "label"}. The default is FALSE, and
#'   the node label would be returned; otherwise, the alias of node label would
#'   be output. The alias of node label is created by adding a prefix
#'   \code{"Node_"} to the node number if the node is an internal node or adding
#'   a prefix \code{"Leaf_"} if the node is a leaf node.
#'
#' @export
#' @return The label of the shared node
#' @author Ruizhu Huang
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
#'           return = "label", use.alias = FALSE)
#'
#' shareNode(node = c('t10','Node_17'), tree = tinyTree,
#'           return = "label", use.alias = FALSE)
#'
#' # -----------------------------------------------------
#' ggtree(tinyTree,branch.length = 'none')+
#' geom_text2(aes(label = node))
#'
#' ## find the node shared by provided node numbers
#' shareNode(node = c(2, 3), tree = tinyTree)

shareNode <- function(tree, node,
                      return = c("number", "label"),
                      use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree is not a phylo object.")
    }

    if (!is.atomic(node)) {
        stop("node is a vector")
    }

    # transfer node label to node number
    if (is.character(node)) {
        node <- transNode(tree, input = node,
                          message = FALSE)
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
    sNode <- as.numeric(df[df$Freq == min(df$Freq), "vec"])

    # final output (node number or label)
    return <- match.arg(return)
    switch(return,
           number = sNode,
           label = transNode(tree = tree, input = sNode,
                             use.alias = use.alias,
                             message = FALSE))

}
