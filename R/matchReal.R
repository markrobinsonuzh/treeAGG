#' match inferred results with the truth
#'
#' \code{matchReal} assigns the inferred nodes (or tips) to the corresponding true
#' signal nodes. For each true signal nodes, the assigned inferred nodes are sorted
#' in increasing order according to the number of their descendant tips
#'
#' @param infN a vector of inferred results (tips & nodes)
#' @param trueN a vector of true results (tips & nodes)
#' @param stree a list of phylo object
#' @param decreasing logical value; indicate the sort order
#'
#' @result a list with each element including the inferred nodes / tips which are
#' assigned to a true node
#'
#' @examples {
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
#' # true nodes
#' tN <- c("Node13","Node2")
#'
#' # inferred nodes & tips
#' iN <- c("Node14","Node16","Node2")
#'
#' mList <- matchReal(infN = iN, trueN = tN, stree = sTree)
#' }
#'

matchReal <- function(infN, trueN, stree, decreasing = FALSE){


 aList <- sapply(trueN, FUN = FindOffspring, stree = stree,
                 simplify = FALSE)
 bList <- sapply(infN, FUN = FindOffspring, stree = stree,
                simplify = FALSE)

 amList <- lapply(aList, FUN = function(x){
    ind <- lapply(bList, FUN = function(y, z){
     ii <- intersect(y, z)
     length(ii) > 0}, z = x)
    sb <- names(bList)[unlist(ind)]})

 bmList <- lapply(amList, sortNode, stree= stree,
                  decreasing= decreasing )
 return(bmList)
}
