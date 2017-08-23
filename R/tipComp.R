
#' Compare the number of descendant tips in a result list to the truth
#'
#' After applying min P-value algorithm (using function find.p.min), a character vector
#' includes the labels of some tips and nodes is suggested to be the location with signals.
#' For each location (node or tip) in the result, we could use function tipComp to compare
#' it with the truth. A pair value (m,n) is returned, where m is the number of descendant
#' tips for the inferred location (tip return 1) and n is the true number of descendant
#' tips for the inferred location.
#'
#' @param rList a list includes all the inferred significant locations
#' (label of tips and nodes)
#' @param trueList a list includes all the true significant locations
#' @return a data frame with 6 columns.
#'
#' Inferred:  the number of tips under the inferred location (tip return 1).
#'
#' Real: the true tip number under a real signal node
#'
#' class:  the real class of inferred location (NoSignal / Cluster number)
#'
#' NumTip:  tip number in the corresponding class
#'
#' Distance: the shortest distance between point (Inferred, Real) to the line defined by function x-y=0; the column "scaleFreq" is derived from
#'
#' scaleFreq: round((Freq*Inferred/NumTip))
#'
#' @export
#'

tipComp<-function(rList,trueList,stree){

  # check object class for each argument
  if (!is.list(rList) & !is.vector(rList)){stop("rList must be a list")}
  if (!is.list(trueList)){stop("trueList must be a list")}
  if (!is.list(stree)){stop("stree must be a list")}

  if(is.list(rList) & length(rList)!=1){
    cat("Error : rList should be a vector or one length List")
  }else{

    # below is to make sure this function could be used in lapply function
    # lapply (eg x <- list(a=1:4,b=2:5); lapply(x,FUN=function(x){
    # print(class(x))}) (Result is : integer instead of list)

    if(length(rList)==1 & is.list(rList)){
      rList <- unlist(rList)
    }
    rList <- as.list(rList)

    ## True set : True signal location (could be nodes, tips, or both)
    ## Inferred set : inferred signal location (could be nodes, tips, or both)

    ## 1) true set on tip level (trueList)
    trueTip <- trueList

    ## 2) inferred set on tip level (infTip)
    # descendant tips for each inferred location (tip return tip as descendant)
    infTip <- lapply(rList,FUN = FindOffspring,stree=Stree)

    ## 3) non-inferred & true tips
    nInfTip <- lapply(trueTip, FUN = function(x){
      setdiff(x,unlist(infTip))
    })

    ind <- unlist(lapply(nInfTip,length))>0
    new_nInfTip <- nInfTip[ind]
    new2_nInfTip <- as.list(unlist(new_nInfTip))
    names(new2_nInfTip) <- NULL

    ## 4) combine inferred & non-inferred set on tip level
    cInfTip <- c(infTip,new2_nInfTip)

    # number of descendant tips for inferred location (inferred column)
    infN <- c(lapply(infTip,FUN = length),
              lapply(new2_nInfTip,FUN = function(x){0}))

    ## 5) match infTip to 1 instead of (trueTip)
    # copy trueTip if there are overlaps
    copyTip <- c(
                # map inferred cluster to branch signal, return number of tips
                # in the signal branch
                 lapply(cInfTip[seq_len(length(infTip))],
                      FUN = function(z){
                         lapply(trueTip,FUN = function(x,y){
                           k<-intersect(x,y)
                           if(length(k)>0){
                             k<-x
                             return(k)}
                           }, y=z)}),
                 # if a signal tip is not inferred, return 1 in real column
                  lapply(cInfTip[-seq_len(length(infTip))],
                        FUN = function(z){
                          lapply(trueTip,FUN = function(x,y){
                            k<-intersect(x,y)
                            if(length(k)>0){
                              k<-y
                              return(k)}
                          }, y=z)})
                 )

    # number of descendant tips with real signal
    LTip0 <- lapply(copyTip,FUN = function(x){lapply(x,FUN = length)})
    LTip1 <- lapply(LTip0,FUN = unlist)
    trueN <- lapply(LTip1,FUN = function(x){
      x <- unlist(x)
      if(all(x == 0)){
        x <- unique(x)
      }else{
        x <- x[x != 0]}})
    nam.trueN <- names(unlist(trueN))
    nam.trueN[nam.trueN == ""] <- "NoSignal"

    # number of tips within different clusters (Nosignal / signal clusters)
    CollectTip <- lapply(Stree,FUN = function(x){x$tip.label})
    AllTip <- unique(unlist(CollectTip))
    nATip <- length(AllTip)

    vTip <- c("NoSignal" = nATip - sum(unlist(lapply(trueTip,length))),
                   unlist(lapply(trueTip,length)))
    NumTip <- vTip[nam.trueN]

    nt <- lapply(trueN,length)

    # construct matrix : 1st colum (number of tips for inferred location)
    #                    2nd colum (number of tips with real signal)
    #                    3rd colum (number of tips within corresponding clusters)

    # calculate overlap frequency

    df0 <- cbind(rep(unlist(infN),nt), unlist(trueN))
    rownames(df0) <- NULL
    cFreq <- FreqMat(cbind(df0,nam.trueN))

    # calculate distance
    Dist <- calcDist(PointX = df0[,1], PointY = df0[,2],
                     LineX = 1, LineY = -1, LineC = 0)

    df <- cbind.data.frame(df0,nam.trueN,cFreq, NumTip, Dist,
                round(cFreq*rep(unlist(infN),nt)/NumTip,digits = 0))
    colnames(df) <- c("Inferred","Real","class","Freq", "NumTip",
                      "Distance","scaleFreq")

    return(df)
  }
}

#' Calculate the frequence of rows with same entries in a matrix
#'
#' FreqMat calculate the frequency of rows with same values in a matrix
#' @param Mat a matrix
#' @return a vector includes the frequency value of each row.
#' @export
#' @examples {
#'  MM <- cbind(c(1,1,1,3,2,2,3),c(0,0,2,1,0,1,1))
#'  kk<-FreqMat(MM)
#'  cbind(MM,kk)
#' }
#'

FreqMat<-function(Mat){
  cTab <- apply(Mat,MARGIN = 1,FUN = function(z){
    apply(Mat,MARGIN = 1,FUN = function(x,y){all(x==y)},y=z)
  })
  cFreq <- apply(cTab,MARGIN = 1,FUN = sum)
  return(cFreq)
}


#' Distance calculation
#'
#' Calculate the shortest distance between a point and a line. The distance equal to
#' 0 if a point is on the line. Points on different sides of the line would have
#' different signs.
#'
#' @param PointX the x value of a point (or points)
#' @param PointX the y value of a point (or points)
#' @param LineX  the parameter of x for line (eg. if line is ax+by+c=0, then LineX = a) )
#' @param LineY  the parameter of y for line (eg. if line is ax+by+c=0, then LineY = b)
#' @param LineC  the parameter of y for line (eg. if line is ax+by+c=0, then LineC = c)


calcDist <- function(PointX,PointY,LineX,LineY,LineC){

  if(length(PointX) != length(PointY)){
    cat("The length of x-axis value is different to that of y-axis for points")
  }
  dist <- LineX*PointX + LineY*PointY + LineC
  return(dist)
}


