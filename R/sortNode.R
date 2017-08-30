#' reorder the nodes according to the number of its desendant tips
#'
#' \code{sortNode} reorder the nodes according to the number of its
#' descendant tips. Here, we let a tip has one descendant tips (itself).
#'
#'
#'@param vNode a vector of tips & nodes
#'@param stree a list of phylo object
#'@param decreasing default true; nodes with more desendant tips come first.
#'
#'@return a vector
#'
#'@examples{
#' library(ape)
#' library(ggtree)

#' n <- 20
#' Lab <- paste("Node",1:(n-1),sep="")
#' tipLab <- paste("T",1:n,sep="")

#' # entire tree
#' set.seed(1)
#' Tree <- rtree(n,tip.label= tipLab)
#' new.Tree <- addNodeLab(treeO = Tree, nodeLab = Lab)
#' sTree <- pruneTree(new.Tree)

#' # PLOT tree
#' ggtree(new.Tree,branch.length = "none")+
#'   geom_tiplab()+
#'   geom_text2(aes(subset=!isTip,label=label),hjust=-0.3)
#'
#' # reorder nodes
#' vn <- c("Node11", "T12", "Node9")
#' sn1 <- sortNode(vn, sTree)
#' sn2 <- sortNode(vn, sTree, decreasing = TRUE)
#'}

sortNode <- function(vNode, stree, decreasing = FALSE){

  if(length(vNode) > 0){
     tList <- sapply(vNode, FindOffspring, stree = stree,
                     simplify = FALSE)
     nNode <- sapply(tList, length)

     # reorder the nodes
     rNode <- sort(nNode, decreasing = decreasing)
     newNode <- names(rNode)
     return(newNode)
   }else{return(vNode)}

}
