#' find the node shared by specified tips / nodes or both
#'
#' \code{shareNode} could find the lowest level node which is shared
#'  by the specified descendants. (the root is the highest level on the tree)
#'
#' @param child a vector of specified tips or/and nodes
#' @param wtree a tree (phylo class)
#' @param stree a list of subtrees which are cut at each node of the wtree
#'
#' @return the label of the shared node
#'
#'
#' @examples {
#' library(ape)
#' n <- 20
#' Lab <- paste("Node",1:(n-1),sep="")
#' tipLab <- paste("T",1:n,sep="")
#'
#' # entire tree
#' Tree <- rtree(n,tip.label= tipLab)
#' Tree <- addNodeLab(treeO = Tree, nodeLab = Lab)
#'
#' # Subtrees
#' sTree <- pruneTree(Tree)
#'
#'
#' # PLOT tree
#' ggtree(Tree,branch.length = "none")+
#' geom_tiplab()+
#' geom_text2(aes(subset=!isTip,label=label),hjust=-0.3)
#'
#' # find share node
#' shareNode(child = c("T11","T4"), wtree =Tree, stree= sTree)
#' shareNode(child = c("T11","T4","Node9"), wtree =Tree, stree= sTree)
#'
#' }

shareNode <- function(child, wtree, stree){

  childList <- lapply(stree, FUN = function(x){c(x$tip.label,
                                                 x$node.label)})
  # find branch which have included the specified child
  branchList <- lapply(childList, FUN = function(x,y){all(y %in% x)},
                       child)
  vbranch <- unlist(branchList)
  vbranch1 <- names(vbranch)[vbranch]

  # select the branch with the lowest number of descendant tips
  Ntip <- lapply(stree[vbranch1], FUN = function(x){x$Nnode+1} )
  Ntip1 <- unlist(Ntip)
  selNode <- names(Ntip1)[Ntip1 == min(Ntip1)]
  return(selNode)
}
