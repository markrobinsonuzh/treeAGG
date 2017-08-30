#' convert inferred results into clusters
#'
#' \code{groupMat} converts the inferred nodes into clusters by combining the tree information.
#' Each inferred node represent a branch. The root node of the branch is the inferred node. All
#' tips in a branch are grouped as a cluster. If the inferred result also include tips, then each
#' inferred tip is a cluster. The tips not included in the inferred branch or inferred tips are
#' grouped as one cluster. Therefore, if the inferred result includes n elements (include both tips
#' and branches), then the final cluster number is n+1.
#'
#' @param vectNT a vector includes inferred nodes (represent branch) and inferred tips
#' @param stree  a list of phylo object
#'
#' @return a vector which assigned a cluster number to each tip of the tree. In total, there are
#' n+1 clusters if the inferred result has length n.
#'
#' @examples {
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
#' infNodes <- c("Node2", "T12")
#'
#' # final result
#' cc <- groupMat(vectNT = infNodes, stree = sTree)
#'
#' }

groupMat <- function(vectNT, stree ){

  if(!is.atomic(vectNT)){
    stop("vectNT need to be an atomic vector")
  }
  if(!(is.list(stree) & class(stree[[1]]) == "phylo")){
    stop("stree need to be a list of phylo object; see function pruneTree
         to create a list of phylo object")
  }

  # inferred clusters
  infList <- sapply(vectNT, FUN = FindOffspring, stree,
                    simplify = FALSE)

  # all tips
  tipList <- lapply(stree, FUN = function(x){x$tip.label})
  allTip <- unique(unlist(tipList))

  # tips no inferred
  noSig <- setdiff(allTip, unlist(infList))

  # list for inferred & non-inferred
  allList <- c(infList, list(noSig))

  # cluster matrix
  clusNum <- as.list(seq_along(allList))
  clusLen <- lapply(allList, length)
  clusList <- mapply(rep, clusNum, clusLen)

  vectClus <- unlist(clusList)
  names(vectClus) <- unlist(allList)

  ff <- vectClus[allTip]
  return(ff)
  }

