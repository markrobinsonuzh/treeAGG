#' Find offsprings
#'
#' Find offsprings for nodes.
#'
#' @param ancestor A charactor vector includes the label of nodes and tips.
#' If tips are included, themselves are returned as offsprings
#' @param stree a list of all the subtrees of a full phylogenetic tree
#' @param only.Tip if TRUE, only descendant tips are returned; otherwise,
#'         both tips and nodes are returned
#' @return a character vector which includes labels of offspring tips
#' @examples \dontrun{
#' library(ape)
#' library(ggtree)
#'
#' n <- 20
#' Lab <- paste("Node",1:(n-1),sep="")
#' tipLab <- paste("T",1:n,sep="")
#'
#' # entire tree
#' set.seed(1)
#' Tree <- rtree(n,tip.label= tipLab)
#' new.Tree <- addNodeLab(treeO = Tree, nodeLab = Lab)
#' sTree <- pruneTree(new.Tree)
#'
#' # PLOT tree
#' ggtree(new.Tree,branch.length = "none")+
#' geom_tiplab()+
#' geom_text2(aes(subset=!isTip,label=label),hjust=-0.3)
#'
#' # find offspring
#' offsp1 <- FindOffspring(ancestor = "Node14", stree = sTree,
#' only.Tip = TRUE)
#' offsp2 <- FindOffspring(ancestor = "Node14", stree = sTree,
#' only.Tip = FALSE)
#' }


elementFind <- function(ancestor,stree, only.Tip = TRUE){

  # check whether all ancestor specified exist in the tree
  UnitAll <- unique(unlist(lapply(stree,
                                  FUN = function(x){c(x$tip.label,
                                                      x$node.label)})))
  notIn <- setdiff(ancestor,UnitAll)

   if(length(notIn) > 0){
    stop(cat(notIn,"cann't be found from the tree"))
   }

  # find offsprings; tips themselves are returned as their offsprings
  node.name <- names(stree)
  nod.res<-ancestor[ancestor %in% node.name]

  if(only.Tip){
    des.res<-unlist(lapply(stree[nod.res]
                           ,FUN = function(x){x$tip.label}))
  }else{
    des.res<-unlist(lapply(stree[nod.res]
                           ,FUN = function(x){c(x$tip.label,x$node.label)}))
  }


  des.final<-c(des.res,setdiff(ancestor,nod.res))
  names(des.final) <- NULL

  return(des.final)
}

FindOffspring <- function(ancestor,stree, only.Tip = TRUE){

  a.List <- as.list(ancestor)
  osList <- lapply(a.List, FUN = elementFind, stree = stree,
                   only.Tip = only.Tip)
  if(length(osList) == 1){
    osList <- unlist(osList)
  }else{ osList <- osList }

  return(osList)

}
