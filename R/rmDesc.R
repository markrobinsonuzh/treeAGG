#' Remove nodes which are the descendants of other nodes
#'
#' \code{rmDesc} removes branches, which are the subbranches of the others.
#'
#' @param node A vector of node labels or node numbers
#' @param tree A list of phylo objects.
#'
#' @export
#'
#' @return A vector of node numbers
#' @examples {
#' data(tinyTree)
#' library(ggtree)
#' # PLOT
#' ggtree(tinyTree,branch.length = 'none') +
#'   geom_text2(aes(label=node), hjust=-.3)
#'
#' # find the shared nodes from the tree plot
#' aV <- c(5, 4, 18, 3, 2)
#' # final result
#' (rn <- rmDesc(node = aV, tree = tinyTree))
#'}
#'

rmDesc <- function(tree, node) {

  if (!inherits(tree, "phylo")) {
    stop("tree: should be a phylo object")
  }

  if (inherits(node, "character")) {
    node <- transNode(tree = tree, input = node)
  } else {
    node <- node
  }

  ListTip <- lapply(as.list(node), FUN = function(x) {
    findOS(ancestor = x, tree = tree, only.Tip = TRUE, self.include = TRUE)
  })
  ListTip1 <- ListTip[!duplicated(ListTip)]
  indList <- lapply(ListTip1, FUN = function(x) {
    aa <- lapply(ListTip1, FUN = function(y, z) {
      all(z %in% y)
    }, z = x)
    atf <- sum(unlist(aa)) == 1
    return(atf)
  })

  ind <- unlist(indList)
  fNT <- node[!duplicated(ListTip)][ind]

  return(fNT)
}
