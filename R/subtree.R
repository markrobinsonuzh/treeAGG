#' prune tree at each node
#'
#' This function prune the tree at each node.
#'
#' @param wtree an object of phylo class
#'
#' @return a list of phylo object
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
#' new.Tree <- addNodeLab(treeO = Tree, nodeLab = Lab)
#' smallTree <- pruneTree(new.Tree)
#' pruneTree(Tree)
#' }

# !!!! modify this to a more stable way
pruneTree <- function(wtree){

  if(is.null(wtree$node.label)){
    stop("node labels of the tree should be specified in advance;
             use addNodeLab to specify the tree nodes")
  }

  small.tree<-ape::subtrees(wtree)
  names(small.tree)<- wtree$node.label

  return(small.tree)
  }
