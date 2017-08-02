#' cut tree into clusters
#'
#' \code{ClusterForm} cut tree into a defined number clusters,
#' and save tips in each cluster as a list
#'
#' @param wtree a phylo tree to be cut;
#' @param nClus number of clusters
#'
#' @return a list with each element representing a cluster.
#' In each element, all tips in that cluster are saved.
#'
#' @example
#' library(ape)
#' n <- 20
#' nodeLab <- paste("Node",1:(n-1),sep="")
#' tipLab <- paste("T",1:n,sep="")
#' # generate an ultrametric tree
#' Tree <- rcoal(n,tip.label= tipLab)
#' Tree$node.label <- nodeLab
#'
#' ClusterForm(Tree,3)

ClusterForm <- function(wtree,nClus){

  # check phylo class
  if(class(wtree) == "phylo"){
    wtree <- ape::as.hclust.phylo(wtree)
  }

  # cut tree into different clusters
  cluster <- cutree(wtree,k=nClus)

  # group tips into clusters and save as a list
  numClus <- as.list(1:nClus)
  names(numClus) <- paste("cluster",1:nClus,sep="")
  ListClus <- lapply(numClus,FUN = function(x){
    names(which(cluster == x))
  })

  return(ListClus)
}


#' Information of clusters
#' \code{ClusterInf} summarize the total number of tips in the cluster
#' and the proportion of counts for each cluster calculated from the real data
#'
#' @param RealDat count table (include at least all tips of the tree)
#' @param tree a phylo tree to be cut;
#' @param nClus number of clusters
#'
#' @return a data frame
#'
#' First column indicates the cluster name
#' Second column indicates the total number of tips in the cluster
#' Third column indicates the proportions of counts for the clusters which are estimated
#' from the real data use \code{dirmult} from package dirmult


ClusterInf <- function(RealDat, wtree, nClus, stree){

  clusList <- ClusterForm(wtree=wtree, nClus)
  # estimate parameters for Dirichlet distribution
  DirMultOutput<-dirmult(data=as.matrix(RealDat))

  ############# tips & clusters proportion ###############
  p.est<-DirMultOutput$pi
  names(p.est) <- names(DirMultOutput$pi)

  # cluster proportion
  p.clus <- lapply(clusList,FUN = function(x){p.est[x]})

  # cluster name & total number of tip & proportion of counts
  ClusInf <- cbind.data.frame("Cluster" = names(clusList),
                              "Freq" = unlist(lapply(clusList,length)),
                              "Proportion" = unlist(lapply(p.clus,sum)))
  return(ClusInf)

}
