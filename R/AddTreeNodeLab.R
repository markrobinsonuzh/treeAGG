#' add node label to the tree
#'
#' This function add labels to all tree nodes.
#'
#' @param treeO an object of phylo class
#' @param nodeLab the label to be used for the tree nodes
#'
#' @return an object of phylo class
#'
#' @examples {
#'
#' library(ape)
#' n <- 20
#' Lab <- paste("Node",1:(n-1),sep="")
#' tipLab <- paste("T",1:n,sep="")
#'
#' # entire tree
#' Tree <- rtree(n,tip.label= tipLab)
#'
#'  # PLOT tree
#' ggtree(Tree,branch.length = "none")+
#' geom_tiplab()+
#' geom_text2(aes(subset=!isTip,label=label),hjust=-0.3)
#'
#' new.Tree <- addNodeLab(treeO = Tree, nodeLab = Lab)
#'
#' # node labels are added
#' ggtree(new.Tree,branch.length = "none")+
#' geom_tiplab()+
#' geom_text2(aes(subset=!isTip,label=label),hjust=-0.3)
#'
#' }

addNodeLab <- function(treeO,nodeLab){

  wtree <- treeO
  wtree$node.label <- nodeLab

  return(wtree)
}




