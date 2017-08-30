#' remove super branches
#'
#' \code{rmSuper} is to remove the branches which have sub-branches or tips.
#' If many branches (represented by its root node) and tips are given as a vector,
#' this function could find and remove the branch which has at least one
#' subbranch or tip existed in the given vector.
#'
#'@param vectNT a vector; a branch is represented by its root node.
#'@param stree a list of phylo object.
#'
#'@examples {
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

#' # find the shared nodes from the tree plot
#' aV <- c("Node5","Node4","Node2","T11","Node18","Node12","T3")

#' # final result
#' rn <- rmSuper(vectNT = aV, stree = sTree)
#'}
#'

rmSuper <- function(vectNT, stree){
  # descendant tips
  tipList <- lapply(as.list(vectNT), FUN = FindOffspring,
                    stree, only.Tip = FALSE)
  names(tipList) <- vectNT

  # number of branches & tips included
  indList <- lapply(tipList, FUN = function(x, y){
     aa <- intersect(x, y)
     la <- length(aa)}, y = vectNT)

  # if has descendant branch or tip (>=2), remove
  ind <- unlist(indList) < 2

  fr <- vectNT[ind]
  return(fr)
}

