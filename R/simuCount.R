#' simulate a count table based on the signals on a tree
#'
#' \code{simuCount} could simulate a count table for all tips and nodes of a tree
#' on two different conditions. The signal on a specified part of a tree could be
#' easily simulated by changing arguments like swapClus or diffClus in \cpde{simuCount}
#' function. The proportions of all tips are estimated from the real data.
#'
#' @param RealDat real data used to estimate the proportion of tips on a tree
#' @param wtree a tree (phylo class) where signal branch could be specified
#' @param nClus the number of clusters which we would like to cut the tree
#'              and finally select some clusters to add signal
#' @param nSam the sample number for each condition
#' @param muNB mean for negative binomial (mu argument in stats::rnbinom)
#'             (to sample library size)
#' @param sizeNB target for number of successful trials (size argument in stats::rnbinom)
#'                (to sample library size)
#' @param swapClus a vector specifies two branches we want to swap their proportions
#'                 (use the root node labels for these two branches); By default,
#'                 NULL is set, which means we don't swap any two branches.
#' @param diffClus a vector specifies branches which we want to add signal
#'                 (differentially abundant) to. The signal is created by directly
#'                 changing their fold change; By default, NULL is set,
#'                 which means we don't change fold change of any branch.
#'                 (use the root node label for each branche)
#'
#' @param FC a numeric vector indicates the fold change. The order of FC and
#' diffClus are matched. If diffClus is set to NULL, FC should be set to NULL too
#'
#'
#' @return a matrix with rows representing tips of the tree and columns representing
#' samples. Samples in condition 1 and condition 2 are named corresponding with
#' the start string "C1_" and string "C2_"
#'
#'
#' @examples {
#' library(dirmult)
#' library("GUniFrac")
#' data(throat.tree)
#' data(throat.otu.tab)
#'
#' Wtree <- addNodeLab(treeO = throat.tree,
#' nodeLab = paste("Node",1:throat.tree$Nnode,sep=""))
#'
#' DF <- simuCount(RealDat = throat.otu.tab,
#' wtree = Wtree,nClus = 40,nSam = 5,
#' muNB = 10000,sizeNB = 5,
#' swapClus = c("cluster1","cluster19"),
#' diffClus = NULL,FC=NULL)
#'}


simuCount <- function(RealDat, wtree, nClus,
                      nSam, muNB, sizeNB,
                      swapClus,diffClus = NULL,FC = NULL){

  # estimate parameters for Dirichlet distribution
  DirMultOutput<- dirmult::dirmult(data=as.matrix(RealDat))

  ############# tips & clusters proportion ###############
  estP<-DirMultOutput$pi
  names(estP) <- names(DirMultOutput$pi)
  # data simulation for condition 1
  Tips <- wtree$tip.label
  p.est <- estP[Tips]   # sort according to the order of tip label

  # cluster name & total number of tip & proportion of counts
  clusList <- ClusterForm(wtree = wtree, nClus = nClus)  # name order 1:nClus
  p.clus <- lapply(clusList,FUN = function(x){p.est[x]})
  ClusInf <- cbind.data.frame("Cluster" = names(clusList),
                              "Freq" = unlist(lapply(clusList,length)),
                              "proportion" = unlist(lapply(p.clus,sum)))


  ############ create count matrix (condition 1) #################

  # parameter alpha for dirichlet distribution
  theta <- DirMultOutput$theta
  gplus <- (1 - theta) / theta
  p.c1 <- p.est
  g.c1 <- p.c1 * gplus

  # ----------- simulate Libary size from NB ----------------
  nSam1 <- nSam[1]
  nSeq <- rnbinom(n=nSam1, mu = muNB, size = sizeNB)

  # ---------------- simulate counts -------------------------
  # the vector of alpha for dirichlet distribution is fixed
  # for each simulation, the proportions are resampled
  # from dirichlet distribution

  # Mp.c1 : sample proportions
  # Mobs.c1 : sample counts
  Mp.c1 <- matrix(0, nSam1, length(g.c1))
  rownames(Mp.c1) <- paste("C1_",1:nSam1,sep="") # tips in columns
  colnames(Mp.c1) <- Tips
  Mobs.c1 <- Mp.c1

  for (i in 1:nSam1) {
    # proportion simulated from Dirichlet distribution
    Mp.c1[i, ] <- dirmult::rdirichlet(n = 1, alpha = g.c1)[1, ]
    # count matrix simulated from multinomial distribution
    Mobs.c1[i, ] <- rmultinom(1, nSeq[i], prob=Mp.c1[i, ])[, 1]
  }

  ############ create count matrix (condition 2) ##########

  # parameters for the condition 2
  p.c2 <- p.est

  # update p.c2 for two different cases : swap or not
  if(!is.null(swapClus)){
    # swap proportions between two clusters
    TipCD1 <- clusList[[swapClus[1]]]
    TipCD2 <- clusList[[swapClus[2]]]
    p.c2[TipCD1] <- (sum(p.est[TipCD2])/sum(p.est[TipCD1]))*p.est[TipCD1]
    p.c2[TipCD2] <- (sum(p.est[TipCD1])/sum(p.est[TipCD2]))*p.est[TipCD2]

     }else{

    # change the proportions of clusters according to the fold change
    # tips selected to do fold change
    clusProp <- lapply(clusList, FUN = function(x){p.est[x]})
    TipFC <- clusProp[diffClus]

    # change the proportion of the selected tips
    FC.sel <- mapply("*",TipFC,FC)
    names(FC.sel) <- NULL
    vFC.sel <- unlist(FC.sel)
    p.c2[names(vFC.sel)] <- vFC.sel

    # normalise proprotions to make sure that they sum up to 1
    p.c2<-p.c2/sum(p.c2)
  }

  # Mp.c2 hold the underlying proportions (condition 2)
  g.c2 <- p.c2 * gplus
  nSam2 <- nSam[2]
  Mp.c2 <- matrix(0, nSam2, length(g.c2))
  rownames(Mp.c2) <- paste("C2_",1:nSam2,sep="")
  colnames(Mp.c2) <- names(g.c2)
  Mobs.c2 <- Mp.c2
  nSeq2 <- rnbinom(n = nSam2, mu = muNB, size = sizeNB)
  for (i in 1:nSam2) {
    Mp.c2[i, ] <- dirmult::rdirichlet(1, g.c2)[1, ]
    Mobs.c2[i, ] <- rmultinom(1, nSeq2[i], prob=Mp.c2[i, ])[, 1]
  }

  # simulated count matrix
  Count<-t(rbind(Mobs.c1,Mobs.c2))


}



#' Aggregate the count of tips to the nodes
#'
#' This function is to calculate the counts at each node. The node count is
#' the sum of the counts of its descendant tips.
#'
#' @param tipTable the tip count table
#' @param wtree a tree (phylo class)
#' @param stree a list of phylo class;
#'              a list of subtrees cut at each node of wtree
#'
#' @return a count table (matrix class) with row representing a tip /node and
#' columns representing samples. Samples are grouped into two conditions.
#' Samples in condition 1 and condition 2 are named corresponding with
#' the start string "C1_" and string "C2_"



nodeCount<-function(tipTable,wtree,stree){
  ## calculate counts for nodes
  nN<-wtree$Nnode
  nNam<-names(stree)

  # calculate counts at nodes
  cNode<-matrix(NA,nrow=nN,ncol=ncol(tipTable))
  for(i in 1:nN){
    node.i<-nNam[i]
    tips.i<-(stree[[node.i]])$tip.label
    cNode[i,]<-apply(tipTable[tips.i,],2,sum)
  }

  rownames(cNode)<-nNam
  return(rbind(tipTable,cNode))
}

