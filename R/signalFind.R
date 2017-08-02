#' Find the signal location (node/tip/both) in a tree
#'
#' Given the known significant differential tips,
#' this functions could find the highest nodes on the tree shared by these tips
#' and aggregate the tips by their shared nodes to make the result short
#'
#' @param TipDiff  the labels of significant differential tips
#' @param wtree    the entire phylogenetic tree which is based to do the aggregation
#' @param stree    a list of subtrees which are cut at the nodes of phylogenetic tree
#' @return a character vector which shows the simplified result : signal nodes and tips
#' which can not be aggregated any more on the tree
#'
#'@example \notrun{
#' library(ape)
#' n <- 20
#' nodeLab <- paste("Node",1:(n-1),sep="")
#' tipLab <- paste("T",1:n,sep="")
#' # entire tree
#' Tree <- rtree(n,tip.label= tipLab)
#' Tree$node.label <- nodeLab
#'
#' # Subtrees
#' sTree <- subtrees(Tree)
#' names(sTree) <- nodeLab
#'
#' # PLOT
#' ggtree(Tree,branch.length = "none") +
#' geom_text2(aes(subset=!isTip, label=label), hjust=-.3) +
#' geom_tiplab()
#'
#' # Differential tips (expect to get "Node4" "Node8" "T1")
#' # find the shared nodes from the tree plot
#' DiffTip <- c("T7","T15","T1","T9","T19")
#'
#' # final result
#' rn <- signalFind(DiffTip, Tree, sTree)
#'}

signalFind <- function(TipDiff,wtree,stree){
  # all tips
  tips<-lapply(stree,FUN = function(x){x$tip.label})

  # find nodes whose descendant tips are in TipDiff
  lt<-lapply(tips,FUN = function(x){all(x %in% TipDiff)})
  ln<-names(which(unlist(lt)==TRUE) )

  # keep signal nodes
  #       1. whose descendant tips are in TipDiff (signal node)
  #       2. Descendant tips between different nodes are exclusive
  Nt<-lapply(stree[ln],FUN = function(z){
         # how many x having tips all included in y
         # return the number of x and named as y
         lo<-lapply(stree[ln],FUN =function(x,y){
                   k <- all(y$tip.label %in% x$tip.label)},z)
         mn<-sum(unlist(lo))
         return(mn)})
  Nt.ul <- unlist(Nt)
  names(Nt.ul) <- ln
  SN <- names(Nt.ul[Nt.ul==1])

  # keep signal nodes and signal tips
  ST <- setdiff(TipDiff, unlist(lapply(stree[SN],
                                 FUN = function(x){x$tip.label})))
  SR <- c(SN,ST)

  return(SR)
}



