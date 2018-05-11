#' remove nodes which are the descendants of other nodes
#'
#' \code{rmMember} is to remove branches, which are the subbranch of the others.
#'
#'@param node a vector of node labels or node numbers
#'@param stree a list of phylo object.
#'
#'@export
#'
#'@return a vector of node numbers
#'@examples {
#' library(ape)
#' n <- 20
#'

#' # entire tree
#' set.seed(1)
#' Tree <- rtree(n)
#'
#' # PLOT
#' ggtree(Tree,branch.length = 'none') +
#'   geom_text2(aes(label=node), hjust=-.3)
#'
#' # find the shared nodes from the tree plot
#' aV <- c(36, 33, 37, 16)

#' # final result
#' rn <- rmMember(node = aV, tree = Tree)
#'}
#'

rmMember <- function(node, tree) {
    
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }
    
    if (inherits(node, "character")) {
        node <- tx_node(tree = tree, input = node)
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
