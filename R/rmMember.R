
#' find the mutually exclusive largest branches
#'
#' \code{rmMember} is to find the largest branches which are mutually exclusive and
#' could cover all given branches. In other words, if we have a lot of branches,
#' in which some branches are the subbranch of the others, we could use this function to
#' remove the subbranches.
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
#'
#' # find the shared nodes from the tree plot
#' aV <- c("Node12","Node14","Node9","T12")

#' # final result
#' rn <- rmMember(vectNT = aV, stree = sTree)
#'}
#'

rmMember <- function(vectNT, stree){

  ListTip <- lapply(as.list(vectNT), FUN = FindOffspring, stree)

  indList <- lapply(ListTip, FUN = function(x){
    aa <- lapply(ListTip, FUN = function(y, z){
      all(z %in% y)}, z = x)
    atf <- sum(unlist(aa)) == 1
    return(atf)
  })

  ind <- unlist(indList)
  fNT <- vectNT[ind]

  return(fNT)
}
