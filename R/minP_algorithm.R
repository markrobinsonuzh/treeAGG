#' Apply minimum P-value algorithm
#'
#' This function applies minimum p-value algorithm to find the branch which
#' is differentially abundant.
#'
#' @param wtree  phylo object; the entire tree for all OTUs
#' @param ResTipNode data frame (include at least :
#'        1. the label of tree nodes and tips as row names
#'        2. adjusted-p value in one column)
#' @param stree a list of all the subtrees of a phylogenetic tree (wtree)
#' @param P.lim the threshold value for adjusted-P
#' @param VarSig column name for testing significance
#' @param VarAGG variable used to do tree aggregation
#' @return a vector of character shows the labels of tips and nodes which
#' have filtered by the minimum P-value algorithm
#' @examples {
#'
#' }
#'

treeAGG<-function(wtree,ResTipNode,stree,P.lim,VarSig,VarAGG){

  # tips & nodes from the whole tree
  all.node<-wtree$node.label
  num.node<-wtree$Nnode
  all.tip <-wtree$tip.label
  num.tip <-length(all.tip)

  # find the location with minimum p-value in the path connecting the root and a tip)
  # tips or nodes with missing adjusted p-value is treated as non-significant
  isSig<-ResTipNode[,VarSig]<= P.lim
  names(isSig)<-rownames(ResTipNode)

  # firstly, set all tips and nodes as FALSE
  keep.tip <-rep(FALSE,num.tip)
  keep.node<-rep(FALSE,num.node)
  names(keep.tip)<-all.tip
  names(keep.node)<-all.node

  # secondly, change significant nodes and tips to TRUE
  keep<-c(keep.tip,keep.node)
  keep[names(isSig)[which(isSig)]]<-TRUE


  for (i in 1:length(all.node)){
    # parent & children
    node.i<-all.node[i]      # parent node
    tree.i<-stree[[node.i]]
    child.i<-setdiff(c(tree.i$node.label
                       ,tree.i$tip.label),
                     node.i)

    # add 1 here to avoid the NA p-value in tips / nodes
    # (NA might be due to the filteration in DESeq or not observed)
    rank <- ResTipNode[,VarAGG]
    names(rank) <- rownames(ResTipNode)
    mRank.child <- min(c(rank[child.i],1),
                    na.rm = TRUE)
    mRank.node <- min(c(rank[node.i],1),
                      na.rm = TRUE)
    if(keep[node.i]){
      # if min value occur at nodes, remove their children
      # otherwise remove the nodes and keep their children
      if(mRank.node>mRank.child){
        keep[node.i]<-FALSE
      }else{
        keep[child.i]<-FALSE}
    }
  }

  keep1<-keep[keep]
  return(names(keep1))
}
