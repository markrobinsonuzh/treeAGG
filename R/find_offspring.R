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
#' rnorm(10)
#' }


FindOffspring<-function(ancestor,stree, only.Tip = TRUE){

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
