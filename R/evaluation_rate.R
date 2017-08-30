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
#' @return a list for TPR, FPR, or FDR (depends on the type selected.)
#'
#'


calcRate<-function(InfMat,TrueLab,type){

  # create a list to store result
  lt <- length(type)
  vl <- vector("list",lt)
  names(vl) <- type
  for(j in 1:lt){
    vl[[j]] <- rep(NA,ncol(InfMat))
  }


  # calculate for each simulation
  for (i in 1:ncol(InfMat)){
    pred.i<-InfMat[,i,drop=FALSE]
    names(pred.i)<-rownames(InfMat)

    # Real positive
    TP.real<-TrueLab[TrueLab==1]
    n.TP.real<-length(TP.real)
    # Real negative
    TN.real<-TrueLab[TrueLab==0]
    n.TN.real<-length(TN.real)

    # reported as positive
    P.test<-pred.i[pred.i==1]
    n.P.test<-length(P.test)

    # reported as positive and is actually positive
    TP.test<-P.test[names(P.test) %in% names(TP.real)]
    n.TP.test<-length(TP.test)
    # reported as positive but is actually negative
    FP.test<-P.test[names(P.test) %in% names(TN.real)]
    n.FP.test<-length(FP.test)


    if("TPR" %in% type){
     # true positive rate
      vl[["TPR"]][i]<-n.TP.test/n.TP.real
    }

    if("FPR" %in% type){
      # false positive rate
      vl[["FPR"]][i]<-n.FP.test/n.TN.real
    }

    if("FDR" %in% type){
     # false discovery rate
      vl[["FDR"]][i]<-n.FP.test/n.P.test
    }

    # Show the process done
    if(!i %% 10){
      cat(i,"out of",ncol(InfMat),"is done \n",sep=" ")
    }
 }
  return(vl)
}


#' Calculate TPR, FPR, FDR after applying min P-value algorithm
#'
#' Calculate TPR (true positive rate), FPR (false positive rate), FDR (false discovery rate)
#' on tip level after min P-value algorithm was applied
#'
#' @param wtree  phylo object; the entire tree for all OTUs
#' @param ResTipNode A list of data frame;  (each data frame includes at least :
#'        1. the label of tree nodes and tips as row names
#'        2. adjusted-p value in one column)#' @param TrueLab numerical vector. Tips as names; value only include 1 (real positive)
#' or 0 (real negative)
#' #' @param stree a list of all the subtrees of a phylogenetic tree (wtree)
#' @param P.lim : the threshold value for adjusted-P
#' @param VarCol: column name for adjusted-P
#' @param TrueLab numerical vector. Tips as names; value only include 1 (real positive)
#' or 0 (real negative)
#' @param type a charactor vector. Select from "TPR" / "FPR" / "FDR" to show which value would be
#' calculated."TPR" (true positive rate), "FPR" (false positive rate), "FDR" (false discovery rate)
#' @return a list of vectors for TPR, FPR, or FDR (depends on the type selected.)
#' @export
#'


calcRateP<-function(wtree,ResTipNode,stree,P.lim,VarCol,TrueLab,type){

  # number of simulation
  nSimRep<-length(ResTipNode)

  # create matrix : row (tips) column (simulation)
  mat1<-matrix(0,nrow = (wtree$Nnode+1),ncol=nSimRep)
  rownames(mat1)<-wtree$tip.label

  for(m in 1:nSimRep){
    ResTipNode.m<-ResTipNode[[m]]
    ResPmin<-find.p.min(wtree,ResTipNode=ResTipNode.m
                     ,stree,P.lim,VarCol)
    Tips<-FindOffspring(ancestor=ResPmin,stree)
    mat1[Tips,m]<-1
    }
    # calculate value listed in type
    PR<-calcRate(InfMat=mat1,TrueLab,type)

   return(PR)
}
