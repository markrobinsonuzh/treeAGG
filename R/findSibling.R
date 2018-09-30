#' find the sibling node
#'
#' \code{findSibling} is to find the sibling node of an \code{input} node.
#'
#' @param tree A phylo object.
#' @param input A numeric or character vector. Node labels or node numbers.
#' @param return "label" or "number". Default is "number", the node number is
#'   returned. If "label", the selected node label is returned.
#' @param use.alias A logical value, TRUE or FALSE. This is an optional argument
#'   that only requried when \code{return = "label"}. The default is FALSE, and
#'   the node label would be returned; otherwise, the alias of node label would
#'   be output. The alias of node label is created by adding a prefix
#'   \code{"Node_"} to the node number if the node is an internal node or adding
#'   a prefix \code{"Leaf_"} if the node is a leaf node.
#'
#' @export
#' @return A vector of node labels or node numbers
#'
#' @examples
#' library(ggtree)
#' data(exTree)
#' ggtree(exTree, branch.length = 'none') %>%
#'     scaleClade(node = 85, scale = 10)+
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#'
#'  findSibling(tree = exTree, input = 86, return = "number")
#'  findSibling(tree = exTree, input = 86, return = "label")
findSibling <- function(tree, input, return = c("number", "label"),
                        use.alias = FALSE){
    # find descendant leaves of input
    inT <- lapply(input,
                  FUN = function(x){
                      findOS(tree = tree,
                             ancestor = x,
                             only.Tip = TRUE)})
    # find the parent node of the input
    pN <- findAncestor(tree = tree, node = input, level = 1)

    # Leaves not included in input
    exT <- lapply(seq_along(pN), FUN = function(x) {
        aT <- findOS(tree = tree, ancestor = pN[x], only.Tip = TRUE)
        setdiff(aT, inT[[x]])
    })

    # replace leaves with their ancestor branch node
    fT <- lapply(exT, FUN = function(x) {
        signalNode(tree = tree, node = x,
                   return = return, use.alias = use.alias)} )
    fT <- unlist(fT)
    # isT <- isTip(tree = tree, input = fT)
    # fn <- fT[!isT]
    return(fT)
}
