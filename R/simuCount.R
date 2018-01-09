#' simulate a count table based on the signals on a tree
#'
#' \code{simuCount} could simulate a count table for all tips and nodes of a tree
#' on two different conditions. The signal on a specified part of a tree could be
#' easily simulated by changing arguments like swapClus or multClus in \cpde{simuCount}
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
#' @param multClus a vector specifies branches which we want to add signal
#'                 (differentially abundant) to. The signal is created by directly
#'                 changing their fold change; By default, NULL is set,
#'                 which means we don't change fold change of any branch.
#'                 (use the root node label for each branche)
#'
#' @param FC.multClus a numeric vector indicates the fold change. The order of FC and
#' multClus are matched. If multClus is set to NULL, FC.multClus should be set to NULL too
#'
#' @param randClus if "random", then the tips are randomly selected from subset ranges.
#' tips are grouped based on the quantiles of the tip proportions. The quantiles are produced
#' correspoinding to the probabilities given in probQ. In each group, a tip is randomly selected.
#' @param probQ numeric vector of probabilities with values in [0, 1]. This probability is used
#' to produce quantile of tip proportions.
#' @param selQ a vector of length 2, which specifys the interval should be used to do sampling.
#'  For example, if selQ = c(1, 2), then the first and second intervals are selected. The interval
#'  are the different ranges created by by the quantiles of tip proportions
#'  
#' @param nsel the number of randomly selected tips, whose proportions are in an interval
#' selected by selQ
#' @param type "swap", "multiple", "random" or NULL. If swap, two selected clusters are
#' selected to swap their proportions. If multiple, multiple clusters are selected to change
#' their proportions according to FC.multClus; if random, 2*nsel tips are randomly selected
#' to change proportions
#' @param numSim a numeric value; the simulation number
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
#'# swap
#' DF <- simuCount(RealDat = throat.otu.tab,
#' wtree = Wtree,nClus = 40,nSam = c(5,5),
#' muNB = 10000,sizeNB = 5,
#' swapClus = c("cluster1","cluster19"),
#' multClus = NULL,FC.multClus=NULL, type = "swap",
#' numSim=1)
#' 
#' # random
#' DF <- simuCount(RealDat = throat.otu.tab,
#' wtree = Wtree,nSam = c(5,5), muNB = 10000,
#' sizeNB = 5, type = "random", probQ = seq(0,1,length.out=6),
#' selQ = c(1, 5), nsel = 2, numSim=1)
#'}


simuCount<- function(RealDat  = NULL, wtree = NULL, nClus = NULL, 
                     nSam = NULL, muNB = NULL, sizeNB = NULL,
                     swapClus = NULL, multClus = NULL, FC.multClus = NULL,
                     type= "swap", probQ= NULL, selQ, nSel, numSim = 1){
  
  
  # estimate parameters for Dirichlet distribution
  DirMultOutput<- dirmult::dirmult(data=as.matrix(RealDat))
  
  ############# tips & clusters proportion ###############
  estP<-DirMultOutput$pi
  names(estP) <- names(DirMultOutput$pi)
  # data simulation for condition 1
  Tips <- wtree$tip.label
  p.est <- estP[Tips]   # sort according to the order of tip label
  
  # parameter alpha for dirichlet distribution
  theta <- DirMultOutput$theta
  gplus <- (1 - theta) / theta
  
  # signal location : clusters or randomly selected tips
  if(type == "swap" | type == "multiple"){
    clusList <- ClusterForm(wtree = wtree, nClus = nClus)  # name order 1:nClus
  }
  
  if(type == "random"){
    ranG <- quantile(p.est, probs = probQ, na.rm = TRUE)
    cutIT <- cut(p.est, breaks = ranG, include.lowest = TRUE)
    p.estG <-  split(x = p.est, f = cutIT)
    # randomly pick one tip from each group
    p1 <- lapply(p.estG[selQ], FUN = function(x){
      sample(x, size = nSel, replace = FALSE)})
  }
  
  
  
  
  
  # cluster name & total number of tip & proportion of counts
  resList <- lapply(seq_len(numSim), FUN = function(j){
    
    ############ create count matrix (condition 1) #################
    Mobs.c1 <- countC1(p.est = p.est, gplus = gplus,
                       nSam = nSam, muNB = muNB, 
                       sizeNB = sizeNB)
    
    ############ create count matrix (condition 2) ##########
    # swap two clusters
    if(type == "swap"){
      
      if(is.null(swapClus)){
        stop("specify two branch to swap via swapClus")}
      
      Mobs.c2 <- countC2_swap2(p.est = p.est, swapClus = swapClus, 
                               clusList = clusList, gplus = gplus,
                               nSam = nSam, muNB = muNB, 
                               sizeNB = sizeNB)
      
    }
    # change fold change of multiple clusters
    if(type == "multiple"){
      if(is.null(multClus)){
        stop("multClus is missing. Please specify branches which would have signals via multClus")
      }
      if(is.null(FC.multClus)){
        stop("FC.multClus is missing. Please specify fold changes for the multiple signal branches")
      }
      Mobs.c2 <- countC2_mult(p.est = p.est, clusList = clusList, 
                              multClus = multClus, FC.multClus = FC.multClus,
                              gplus = gplus, nSam = nSam, muNB = muNB, 
                              sizeNB = sizeNB)
    }
    
    # change fold change in tips randomly selected
    if(type == "random"){
      if(is.null(probQ)){
        stop("probQ is missing.")
      }
      Mobs.c2 <- countC2_random(p.est = p.est, p1 = p1, 
                                gplus = gplus, nSam = nSam, 
                                muNB = muNB, sizeNB = sizeNB)
      
    }
    
    # no fold change in the whole tree
    if(is.null(type)){
      Mobs.c2 <- countC1(p.est = p.est, gplus = gplus,
                         nSam = nSam, muNB = muNB, 
                         sizeNB = sizeNB)
    }
    
    
    # simulated count matrix
    mat<-list("count"=t(rbind(Mobs.c1,Mobs.c2[[1]])),
              "signal"= Mobs.c2[[2]])
    return(mat)
  })
  
  tabCount <- lapply(resList,FUN = function(x){x$count})
  sigList <- lapply(resList, FUN = function(x){x$signal})
  sigMat <- sigList[!duplicated(sigList)]
  
  if(length(tabCount) == 1){
    tabCount <- do.call(rbind, tabCount)
  }else{ tabCount <- tabCount }
  
  if(length(sigMat) == 1){
    sigMat <- do.call(rbind, sigMat)
  }else{ sigMat <- sigMat }
  
  final <- list(table = tabCount, signal = sigMat)
  return(final)
}




countC1 <-  function(p.est, gplus, nSam, muNB, sizeNB){
  
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
  colnames(Mp.c1) <- names(p.est)
  Mobs.c1 <- Mp.c1
  
  for (i in 1:nSam1) {
    # proportion simulated from Dirichlet distribution
    Mp.c1[i, ] <- dirmult::rdirichlet(n = 1, alpha = g.c1)[1, ]
    # count matrix simulated from multinomial distribution
    Mobs.c1[i, ] <- rmultinom(1, nSeq[i], prob=Mp.c1[i, ])[, 1]
  }
  
  return(Mobs.c1)
}

countC2_swap2 <- function(p.est, swapClus, clusList, gplus,
                          nSam, muNB, sizeNB){
  
  ############ create count matrix (condition 2) ##########
  # parameters for the condition 2
  p.c2 <- p.est
  
  # swap proportions between two clusters
  TipCD1 <- clusList[[swapClus[1]]]
  TipCD2 <- clusList[[swapClus[2]]]
  p.c2[TipCD1] <- (sum(p.est[TipCD2])/sum(p.est[TipCD1]))*p.est[TipCD1]
  p.c2[TipCD2] <- (sum(p.est[TipCD1])/sum(p.est[TipCD2]))*p.est[TipCD2]
  
  
  # Mp.c2 hold the underlying proportions (condition 2)
  g.c2 <- p.c2 * gplus
  nSam2 <- nSam[2]
  Mp.c2 <- matrix(0, nSam2, length(g.c2))
  rownames(Mp.c2) <- paste("C2_",1:nSam2,sep="")
  colnames(Mp.c2) <- names(p.est)
  Mobs.c2 <- Mp.c2
  nSeq2 <- rnbinom(n = nSam2, mu = muNB, size = sizeNB)
  for (i in 1:nSam2) {
    Mp.c2[i, ] <- dirmult::rdirichlet(1, g.c2)[1, ]
    Mobs.c2[i, ] <- rmultinom(1, nSeq2[i], prob=Mp.c2[i, ])[, 1]
  }
  signalMat <- matrix(data = c(p.est[unlist(clusList[swapClus],use.names = FALSE)],
                               p.c2[unlist(clusList[swapClus],use.names = FALSE)]),
                      nrow = 2, byrow = TRUE, 
                      dimnames = list("condition"=c("C1", "C2"), 
                                      "tips"=c(unlist(clusList[swapClus],
                                                      use.names = FALSE))) )
  Mat <- list("Mobs.c2" = Mobs.c2, "signal" = signalMat)
  return(Mat)
}


countC2_mult <- function(p.est, clusList, multClus, FC.multClus,
                         gplus, nSam, muNB, sizeNB){
  
  ############ create count matrix (condition 2) ##########
  # parameters for the condition 2
  p.c2 <- p.est
  
  clusProp <- lapply(clusList, FUN = function(x){p.est[x]})
  TipFC <- clusProp[multClus]
  
  # change the proportion of the selected tips
  FC.sel <- mapply("*",TipFC,FC.multClus)
  names(FC.sel) <- NULL
  vFC.sel <- unlist(FC.sel)
  p.c2[names(vFC.sel)] <- vFC.sel
  
  
  # Mp.c2 hold the underlying proportions (condition 2)
  g.c2 <- p.c2 * gplus
  nSam2 <- nSam[2]
  Mp.c2 <- matrix(0, nSam2, length(g.c2))
  rownames(Mp.c2) <- paste("C2_",1:nSam2,sep="")
  colnames(Mp.c2) <- names(p.est)
  Mobs.c2 <- Mp.c2
  nSeq2 <- rnbinom(n = nSam2, mu = muNB, size = sizeNB)
  for (i in 1:nSam2) {
    Mp.c2[i, ] <- dirmult::rdirichlet(1, g.c2)[1, ]
    Mobs.c2[i, ] <- rmultinom(1, nSeq2[i], prob=Mp.c2[i, ])[, 1]
  }
  
  
  signalMat <- matrix(data = c(p.est[unlist(clusList[multClus],use.names = FALSE)],
                               p.c2[unlist(clusList[multClus],use.names = FALSE)]),
                      nrow = 2, byrow = TRUE, 
                      dimnames = list("condition"=c("C1", "C2"), 
                                      "tips"=c(unlist(clusList[multClus],
                                                      use.names = FALSE))) )
  Mat <- list("Mobs.c2" = Mobs.c2, "signal" = signalMat)
  return(Mat)
}

countC2_random<- function(p.est, p1, gplus, nSam,
                          muNB, sizeNB){
  # calculate fold change
  p.c2 <- p.est
  p.c2[names(p1[[1]])] <- p.est[names(p1[[2]])]
  p.c2[names(p1[[2]])] <- p.est[names(p1[[1]])]
  
  # Mp.c2 hold the underlying proportions (condition 2)
  g.c2 <- p.c2 * gplus
  nSam2 <- nSam[2]
  Mp.c2 <- matrix(0, nSam2, length(g.c2))
  rownames(Mp.c2) <- paste("C2_",1:nSam2,sep="")
  colnames(Mp.c2) <- names(p.est)
  Mobs.c2 <- Mp.c2
  nSeq2 <- rnbinom(n = nSam2, mu = muNB, size = sizeNB)
  for (i in 1:nSam2) {
    Mp.c2[i, ] <- dirmult::rdirichlet(1, g.c2)[1, ]
    Mobs.c2[i, ] <- rmultinom(1, nSeq2[i], prob=Mp.c2[i, ])[, 1]
  }
  
  
  signalMat <- matrix(data = c(p.est[c(names(p1[[1]]),names(p1[[2]]))],
                               p.c2[c(names(p1[[1]]),names(p1[[2]]))]),
                      nrow = 2, byrow = TRUE, 
                      dimnames = list("condition"=c("C1", "C2"), 
                                      "tips"=c(names(p1[[1]]),names(p1[[2]]))))
  
  Mat <- list("Mobs.c2" = Mobs.c2, "signal" = signalMat)
  return(Mat)
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
  nNam<-wtree$node.label
  
  # calculate counts at nodes
  cNode<-matrix(NA,nrow=nN,ncol=ncol(tipTable))
  rownames(cNode)<-nNam
  for(i in 1:nN){
    node.i<-nNam[i]
    tips.i<-(stree[[node.i]])$tip.label
    cNode[i,]<-apply(tipTable[tips.i,],2,sum)
  }
  colnames(cNode) <- colnames(tipTable)
  return(rbind(tipTable,cNode))
}




