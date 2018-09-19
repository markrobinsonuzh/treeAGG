#' find the sibling node
#'
#' \code{findSibling} is to find the sibling node of an \code{input} node.
#'
#' @param tree A phylo object.
#' @param input A numeric or character vector. Node labels or node numbers.
#' @param label A logical value, TRUE or FALSE. If TRUE, node labels are
#' returned.
#'
#' @export
#' @return A vector of node labels or node numbers
#'
#' @examples
#' library(ggtree)
#' data(exTree)
#' ggtree(exTree, branch.length = 'none') +
#'  geom_text2(aes(label = label))
#'
#'  findSibling(tree = exTree, input = 75)
findSibling <- function(tree, input, label = FALSE){
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
        signalNode(tree = tree, node = x, label = FALSE)} )
    fT <- unlist(fT)
    # isT <- isTip(tree = tree, input = fT)
    # fn <- fT[!isT]
    return(fT)
}
