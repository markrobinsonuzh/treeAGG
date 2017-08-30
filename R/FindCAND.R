#' list all possible best results
#' \code{FindCAND} is to find results which are close to the truth
#'
#'@param iNode a vector of inferred nodes
#'@param stree a list of phylo object
#'@param trueN a vector of true nodes (tips)
#'
#'@return a list which lists all found results. Each element includes a vector
#' of kept inferred nodes (tips)
#'
#'@examples{
#' library(ape)
#' n <- 20
#' Lab <- paste("Node",1:(n-1),sep="")
#' tipLab <- paste("T",1:n,sep="")

#' # entire tree
#' set.seed(1)
#' Tree <- rtree(n,tip.label= tipLab)
#' Tree <- addNodeLab(treeO = Tree, nodeLab = Lab)

#' # Subtrees
#' sTree <- pruneTree(Tree)


#' # PLOT
#' ggtree(Tree,branch.length = "none") +
#'   geom_text2(aes(subset=!isTip, label=label), hjust=-.3) +
#'   geom_tiplab()
#'
#' # inferred nodes & tips
#' infNodes <- c("Node2", "Node3", "Node4", "Node7",
#'               "Node8", "Node9", "Node10", "Node11", "Node12",
#'               "Node18", "Node15", "T12")
#' trueNode <- c("Node9","Node17")

#' # final result
#' rn <- FindCAND(iNode = infNodes, stree = sTree, trueN)
#'}
#'

FindCAND <- function(iNode, stree, trueN){

  # signal tips
  stList <- lapply(as.list(trueN), FUN = FindOffspring, stree)
  names(stList) <- trueN

  # non-signal tips
  atList <- lapply(stree, FUN = function(x){x$tip.label})
  nt <- setdiff(unlist(atList), unlist(stList))

  # inferred results
  resList <- sapply(iNode, FUN = FindOffspring, stree)

  # For each signal branch, find the corresponding inferred results
  # which include this signal branch.
  siG <- lapply(stList, FUN = function(x){
          indList <- lapply(resList, FUN = function(y, z){
          isc <- intersect(y,z)
          length(isc) > 0 }, z = x)
          ind <- unlist(indList)
          names(resList)[ind]})

  # the best clustering result below signal node (target node)
  stnList <- lapply(as.list(trueN), FindOffspring, stree,
                      only.Tip = FALSE)
  names(stnList) <- trueN
  # inferred result & significant & below the target node
  lowList <- lapply(stnList, FUN = function(x, y){ y[y %in% x] },
                    y = iNode) # result with signal & below target nodes
  # exclusive & below target node & closest to the target node
  udList <- lapply(lowList, FUN = rmMember, stree)

  # start removing inferred node from the top until it (udList)

  # inferred result & significant & above the target node
  exNode <- setdiff(unlist(siG),unlist(lowList)) # this order is perfect
  # inferred result & significant
  siGNode <- unlist(siG)
  # remove one by one
  if(length(exNode) > 1){
    rmList <- sapply(seq_along(exNode), FUN = function(x){exNode[1:x]})
  }else{rmList <- list(NULL)}

  aftList <- lapply(rmList, FUN = function(x, y){setdiff(y, x)}, y = siGNode)



  # non-signal area : take the node in the lowest tree level and remove
  # their super branch (this is to lower the false positive)
  nsiGNode <- setdiff(iNode, siGNode)
  fNSG <- rmSuper(vectNT = nsiGNode, stree )

  # signal area : select the inferred result which is more similar to the truth
  #
  #   1. above : remove inferred nodes one by one to move closer from top to the true signal node
  #   2. below : keep the inferred nodes which are closest to the true signal (mutually exclusive)
  #   3. if an inferred but non-signal branch is included in an inferred node which covers the
  #      signal area, the inferred but non-signal branch is integrated to the inferred node with signal.
  #      (exclusive requirement)
  aftExList <- lapply(aftList, FUN = function(x, y){
                 u <- c(x, y)
                f <- rmMember(vectNT = u, stree)
                  }, y = fNSG)

  return(aftExList)
  }



