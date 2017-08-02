#' convert tree structure to matrix structure
#'
#' \code{treeLevel} convert a tree structure into a matrix structure with
#' tips in the rows and tree levels in columns.
#'
#'@param wtree the entire tree (phylo class)
#'@param stree a list of subtrees, each of which is cut at a node of wtree
#'
#'@return a matrix. each tip occupy a row.
#'
#' each row includes all nodes in the path connecting the corresponding tip
#' and the root. Columns from left to right includes tree level from high to low.
#' The first column the root label.
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
#' Mtree <- treeLevel(wtree = Tree, stree = sTree)
#'}





treeLevel <- function(wtree,stree){

  # all tips in the whole tree
  tips <- wtree$tip.label
  L.tip <- as.list(tips)
  names(L.tip) <- tips

  # for each OTU : find all nodes in the path connecting the OTU and the root
  #                (included the root )
  L.node <- lapply(L.tip,FUN = function(x){
    nn <-lapply(stree,FUN = function(y,z){
      ind <- z %in% y$tip.label
     },z=x)
    un <- unlist(nn)
    unt <- names(un)[un]
    return(unt)
    })

  # reorder each list element according to the number of descendant tips (decreasing order)
  # nodes with more descendant tips come first
  L1.node <- lapply(L.node, FUN = function(x){
    ln <- lapply(stree[x],FUN = function(y){length(y$tip.label)})
    uln <- unlist(ln)
    oln <- names(uln)[order(uln,decreasing = TRUE)]
  })

  # change from list structure to matrix structure
  # NA are added to the end of elements with short length
  max.length <- max(sapply(L1.node,length))
  L2.node <- lapply(L1.node, FUN = function(x){c(x,rep(NA,max.length - length(x)))})
  M.node <- do.call("rbind",L2.node)
  rownames(M.node) <- names(L1.node)
  colnames(M.node) <- paste("Rank",1:ncol(M.node),sep="")

  return(M.node)
}


