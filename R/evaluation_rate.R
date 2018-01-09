#' Calculate TPR, FPR, FDR on tips level
#'
#' Calculate TPR (true positive rate), FPR (false positive rate), FDR (false discovery rate)
#'
#' @param InfMat a matrix (rows: tree tips; columns: simulations.)
#' The entry value only includes 1 (tested to be positive) or 0 (tested to be negative).
#' @param TrueLab a numerical vector. Tips as names; value only include 1 (real positive)
#' or 0 (real negative)
#' @param type a charactor vector. Select from "TPR" / "FPR" / "FDR" to show which value would be
#' calculated."TPR" (true positive rate), "FPR" (false positive rate), "FDR" (false discovery rate)
#' @param level a numeric value as the threshold value of adjusted p-value. Results with value
#' above this level are categorized as non-significant
#' @return  a matrix with columns for TPR, FPR, or FDR (depends on the type selected.) and rows
#' for simulations
#'
#' @examples {
#' library(dirmult)
#' library("GUniFrac")
#' data(throat.tree)
#' data(throat.otu.tab)
#'
#' # tree
#' Lab <- paste("Node",1:throat.tree$Nnode,sep="")
#' Wtree <- addNodeLab(treeO = throat.tree, nodeLab = Lab)
#' Stree <- pruneTree(wtree = Wtree)
#' clusList <- ClusterForm(wtree = Wtree, nClus = 60)
#'
#' # count table for tips
#' DF.tip <- simuCount(RealDat = throat.otu.tab,
#' wtree = Wtree,nClus = 40,nSam = 5,
#' muNB = 10000,sizeNB = 5,
#' swapClus = c("cluster1","cluster19"),
#' diffClus = NULL,FC=NULL)
#'
#' # count table for tips and nodes
#' DF <- nodeCount(tipTable = DF.tip, wtree = Wtree, stree = Stree)
#'
#' # differential analysis
#' isNode <- substring(rownames(DF), 1, 4) == "Node"
#' resDF <- Redge(countTab = DF, nSam = c(5,5), isTip = !isNode,
#' isAnalyze = rep(TRUE, nrow(DF)), prior.count =5, normalize=TRUE)
#'
#' # true class label
#' sigTip <- unlist(clusList[swapClus])
#' allTip <- Wtree$tip.label
#' Tlabel <- rep(0, length(allTip))
#' names(Tlabel) <- allTip
#' Tlabel[sigTip] <- 1
#'
#' # one simulation
#'  Rate1 <- calcRateP(wtree = throatTree, resModel = mod.edgeR,
#'  stree = smallTree, P.lim = x, VarSig = "FDR",
#'  VarAGG = "FDR", TrueLab=Tlabel, type = c("FDR","TPR"))})
#'
#'  # multiple simulation
#'
#'  }


calcRate<-function(InfMat, TrueLab, type, level, mean){


  # Real positive
  TP.real<-TrueLab[TrueLab==1]
  n.TP.real<-length(TP.real)
  # Real negative
  TN.real<-TrueLab[TrueLab==0]
  n.TN.real<-length(TN.real)

  # calculate for each simulation
  lvv <- lapply(seq_len(ncol(InfMat)), FUN = function(x){
    pred.i<-InfMat[,x,drop=FALSE]
    names(pred.i)<-rownames(InfMat)

    # reported as positive
    P.test<-pred.i[pred.i==1]
    n.P.test<-length(P.test)

    # reported as positive and is actually positive
    TP.test<-P.test[names(P.test) %in% names(TP.real)]
    n.TP.test<-length(TP.test)
    # reported as positive but is actually negative
    FP.test<-P.test[names(P.test) %in% names(TN.real)]
    n.FP.test<-length(FP.test)

    vv <- c("TPR" = n.TP.test/n.TP.real,
            "FPR" = n.FP.test/n.TN.real,
            "FDR" = n.FP.test/n.P.test,
            "level" = level)
    sv <- vv[c(type, "level")]
    return(sv)
  })

  rvv <- do.call(rbind, lvv)
  mvv <- c(apply(rvv[, type, drop=FALSE], 2, mean, na.rm = TRUE),
           "level" = level)
  if(mean){
    return(mvv)
    }else{
    return(rvv)
  }

}




#' Calculate TPR, FPR, FDR after applying min P-value algorithm
#'
#' Calculate TPR (true positive rate), FPR (false positive rate), FDR (false discovery rate)
#' on tip level after min P-value algorithm was applied
#'
#' @param wtree  phylo object; the entire tree for all OTUs
#' @param resModel A list of data frame;  (each data frame includes at least :
#'        1. the label of tree nodes and tips as row names
#'        2. adjusted-p value in one column)#' @param TrueLab numerical vector. Tips as names; value only include 1 (real positive)
#' or 0 (real negative)
#' @param stree a list of all the subtrees of a phylogenetic tree (wtree)
#' @param P.lim : the threshold value for adjusted-P
#' @param TrueLab numerical vector. Tips as names; value only include 1 (real positive)
#' or 0 (real negative)
#' @param type a charactor vector. Select from "TPR" / "FPR" / "FDR" to show which value would be
#' calculated."TPR" (true positive rate), "FPR" (false positive rate), "FDR" (false discovery rate)
#' @param VarSig a character specifys the name of the column which is used to do statistical
#' significance test.
#' @param VarAGG na character specifys the name of the column which is used to do tree aggregation
#' @param mean  a logical value. If TRUE, then the mean of "TPR" / "FPR" / "FDR" on the same FDR
#' level are returned.
#'  @return a list of vectors for TPR, FPR, or FDR (depends on the type selected.)
#'
#'


calcRateP<-function(wtree, stree, resModel, P.lim,
                    TrueLab, type, VarSig, VarAGG,
                    mean = FALSE){

  # number of simulation
  nSimRep <-length(resModel)

  # create matrix : row (tips) column (simulation)
  mat1<-matrix(0,nrow = length(wtree$tip.label),
               ncol=nSimRep)
  rownames(mat1)<-wtree$tip.label

  lat <- lapply(seq_len(nSimRep), FUN = function(x){
    if(is.null(VarSig) && is.null(VarAGG)){
      ResPmin <- resModel[[x]]
    }else{
      resModel.x<-resModel[[x]]
      ResPmin<-treeAGG(wtree,ResTipNode = resModel.x,
                       stree, P.lim, VarSig, VarAGG )
    }

    Tips<-FindOffspring(ancestor=ResPmin,stree)
    m1 <- mat1[, x, drop = FALSE]
    m1[unlist(Tips),] <- 1
    return(m1)
  })

  rat <- do.call(cbind, lat)
  PR<-calcRate(InfMat=rat,TrueLab,type, level = P.lim, mean)

   return(PR)
}
