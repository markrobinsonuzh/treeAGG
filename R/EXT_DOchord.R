

#' produce a matrix with the rows as the inferred result
#' and the columns as the true result.
#' each entry of the matrix is the size of the intersection
#' between inferred result and the true result


choreData <- function(infN, trueN, stree, decreasing){
  # true signal location
  trueN1 <- sortNode(trueN, stree , decreasing)
  tTip <- sapply(trueN1, FindOffspring, stree = stree,
                 simplify = FALSE)
  aTip <- unlist(lapply(stree, FUN = function(x){x$tip.label}))
  oTip  <- list(setdiff(unique(aTip), unlist(tTip)))
  realClass <- c(tTip, NoSignal = oTip)

  # inferred results
  # inferred -- signal
  mList <- matchReal(infN, trueN1 , stree, decreasing)
  # inferred - no signal
  ifN.n<- setdiff(infN, unlist(mList))

  if(length(ifN.n) > 0){
    ifN.n1 <- sortNode(vNode = ifN.n, stree, decreasing)
    }else{
    ifN.n1 <- ifN.n
    }
    #names(ifN.n1) <- rep("falsePos", length(ifN.n1))
    iVect <- c(unlist(mList, use.names = FALSE), ifN.n1)
    # if more than one signal branches are included in one inferred signal set
    iVect <- unique(iVect)
    infTip  <- sapply(iVect, FindOffspring, stree = stree,
                      simplify = FALSE)
    infNeg <- list("NoSignal" = setdiff(unique(aTip), unlist(infTip)))
    infR <- c(infTip, infNeg)



  xx <- lapply(realClass, FUN = function(x){
        si <- lapply(infR, FUN = function(y, z){
        ii <- intersect(y, z)
        li <- length(ii)
        return(li)}, z = x)

        ui <- unlist(si)
        return(ui) })

  rbx <- do.call(rbind, xx)
  Mat <- t(rbx)

  return(Mat)
}

# calculate the location of each polygon item
axisValue <- function(df, rIn, rOut, itv, spPCT,
                      gap.a, gap.b, ds, left = TRUE){
  # row / colum sum
  rNum <- apply(df, 1, sum)
  cNum <- apply(df, 2, sum)

  # row -- left  |  column -- right
  if(left){
    Num <- rNum
  }else{
    Num <- cNum
    }

  # leave gap between left/right part and the middle axis
  ppi <- pi - (gap.a+gap.b)

  # leave space between different sets
  if(length(Num)>1){
    sp <- ppi * spPCT / (length(Num)-1)
  }else{
    spPCT <- 0
    sp <- 0}
  space <- rep(sp, length(Num)-1)


  # polygon size for each set depends on the set proportion
  theta <- cumsum(Num)/sum(Num) * ppi * (1-spPCT)
  theta1 <- c(0, head(theta, -1) + cumsum(space))+gap.a
  theta2 <- c(theta[1], tail(theta, -1) + cumsum(space))+gap.a

  # signal -- differential abundant
  vSig <- names(Num) != "NoSignal"

  xx <- c()
  yy <- c()
  set <- c()
  isSig <- c()

  for (i in seq_len(length(Num))){

    if(!vSig[i]){
      rIn <- rIn - ds
      rOut <- rOut - ds}

    # minor theta
      ai <- itv[itv >= theta1[i] & itv <= theta2[i]]
      # start, minor, end theta
      bi <- c(theta1[i], ai, theta2[i])
      ci <- sort(unique(bi), decreasing = FALSE)

      xx <- c(xx, rIn * sin(ci), rOut * sin(rev(ci)),
              rIn * sin(ci[1]))
      yy <- c(yy, rIn * cos(ci), rOut * cos(rev(ci)),
              rIn * cos(ci[1]))
      set <- c(set, rep(i, length(ci) * 2 + 1))

      isSig <- c(isSig, rep(vSig[i], length(ci) * 2 + 1))
  }

  df <- data.frame(x = xx, y= yy, set = set, isSig = isSig,
                   stringsAsFactors = FALSE)
  return(df)
}

# calculate the end location of the link (end at the left)

linkEnd <- function(df, spPCT.l, gap.a, gap.b){

  # row sum (plot : left part of link)
  Num <- apply(df, 1, sum)
  vSig <- names(Num) != "NoSignal"

  # gap
  ppi <- pi - (gap.a + gap.b)

  # space
  if(nrow(df)==1){
    spPCT.l <- 0
    sp <- 0
    space <- rep(sp, nrow(df) - 1)
    }else{
      sp <- ppi * spPCT.l / (nrow(df) -1)
      space <- rep(sp, nrow(df) - 1)
    }


    # angle limit for each set
  theta <- cumsum(Num)/sum(Num) * ppi * (1-spPCT.l)
  theta.1 <- c(0, head(theta, -1) + cumsum(space)) + gap.a # lower limit
  theta.2 <- c(theta[1], tail(theta, -1) + cumsum(space))+gap.a # upper limit

  theta.e <- c()
  true.e <- c()
  est.e <- c()
  isSig <- c()
  lab.e <- c()
  for(i in seq_len(nrow(df))){
    ri <- df[i, ]
    j <- which(ri != 0)  # corresponding true group
    ri2 <- ri[j]
    pri <- cumsum(ri2)/sum(ri2)

    # angle limits (correspond to different true groups
    # within an estimated set)
    wi <- theta.2[i] - theta.1[i]
    ti <- theta.1[i] + wi * pri
    ti2 <- rep(ti, each = 2)
    nti <- c(theta.1[i], head(ti2,-1))  # include start point
    theta.e <- c(theta.e, nti)

    # indicator : whether the set is inferred as differential abundant
    sig.i <- vSig[i]
    isSig <- c(isSig, rep(sig.i, length(nti)))

    # true group
    j2 <- rep(j, each=2)
    #nj <- c(j[1], head(j2, -1)) # the start point has the same class with the second point
    true.e <- c(true.e, j2)

    # estimated set
    est.e <- c(est.e, rep(i, length(j2)))
    lab.e <- c(lab.e, rep(names(Num)[i], length(j2)))
  }


  df.e <- data.frame( theta.end = theta.e, trueGrp = true.e,
                      estGrp = est.e, isSig.end = isSig,
                      lab.end = lab.e, stringsAsFactors = FALSE)
  df.e1 <- df.e[with(df.e, order(trueGrp, theta.end)),]

  return(df.e1)
}


# calculate the start location of the each link (start from the right)
linkStart <- function(df, spPCT.r, gap.a, gap.b){

  # column sum
  Num <- apply(df, 2, sum)
  vSig <- names(Num) != "NoSignal"

  # gap between left/right part and middle axis
  ppi <- pi - (gap.a + gap.b)

  # space between polygons for different groups
  if(ncol(df)==1){
    spPCT.r <- 0
    sp <- 0
    space <- rep(sp, ncol(df) - 1)
  }else{
    sp <- ppi * spPCT.r / (ncol(df) -1)
    space <- rep(sp, ncol(df) - 1)
  }
  # angle limit for each class
  theta <- cumsum(Num)/sum(Num) * ppi * (1-spPCT.r)
  theta.1 <- c(0, head(theta, -1) + cumsum(space)) + gap.a
  theta.2 <- c(theta[1], tail(theta, -1) + cumsum(space)) + gap.a

  theta.s <- c()
  isSig <- c()
  trG <- c()
  lab.s <- c()

  for(i in seq_len(ncol(df))){
    ri <- df[,i ]
    j <- which(ri != 0)  # corresponding true class
    ri2 <- ri[j]
    pri <- cumsum(ri2)/sum(ri2)

    # angle limits for different groups
    wi <- theta.2[i] - theta.1[i]
    ti <- theta.1[i] + wi * pri
    ti2 <- rep(ti, each = 2)
    nti <- c(theta.1[i], head(ti2, -1))  # include start point
    theta.s <- c(theta.s, nti)

    # indicator : whether the group is differential abundant
    sig.i <- vSig[i]
    isSig <- c(isSig, rep(sig.i, length(nti)))

    # true group
    trG <- c(trG, rep(i, length(nti)))
    lab.s <- c(lab.s, rep(names(Num)[i], length(nti)))
  }

  df <- data.frame(theta.start = theta.s, isSig.start = isSig,
                   trueGrp = trG, lab.start = lab.s,
                   stringsAsFactors = FALSE)
  return(df)
}

# to create link between real and estimated result
link <- function(df, rR, rL, ds ){


  sq <- seq(0, 180, length.out = 361) * (pi/180)
  k <- seq(0, 1, by = 0.02)

  xlink <- c()
  ylink <- c()
  set <- c()
  group <- c()

  for(i in seq_len(nrow(df)-1)){
    # group
    tg1 <- df$trueGrp[i]
    tg2 <- df$trueGrp[i+1]
    eg1 <- df$estGrp[i]
    eg2 <- df$estGrp[i+1]

    # theta
    te1 <- df$theta.end[i]
    te2 <- df$theta.end[i+1]
    ts1 <- df$theta.start[i]
    ts2 <- df$theta.start[i+1]

    # significant indicator
    sig.ei <- df$isSig.end[i]
    sig.si <- df$isSig.start[i]

    if(!sig.ei){
      rL.i <- rL - ds
    }else{rL.i <- rL}
    if(!sig.si){
      rR.i <- rR - ds
    }else{rR.i <- rR}

    if(tg1 == tg2 & eg1 == eg2){

      # if non- differential, short radius
      # end (left)
      sqe <- c(te1, sq[sq <= te2 & sq >= te1], te2)
      xe <- - rL.i * sin(sqe)
      ye <- rL.i * cos(sqe)
      xe1 <- - rL.i * sin(sqe[1])
      ye1 <- rL.i * cos(sqe[1])
      xe2 <- - rL.i * sin(tail(sqe, 1))
      ye2 <- rL.i * cos(tail(sqe, 1))


      # start (right)
      sqs <- c(ts1, sq[sq <= ts2 & sq >= ts1], ts2)
      xs <- rR.i * sin(sqs)
      ys <- rR.i * cos(sqs)
      xs1 <-  rR.i * sin(sqs[1])
      ys1 <- rR.i * cos(sqs[1])
      xs2 <-  rR.i * sin( tail(sqs,1) )
      ys2 <- rR.i * cos( tail(sqs,1) )

      # above
      xa <- c(xe1 * (1-k)^2 + xs1 * k^2)
      ya <- c(ye1 * (1-k)^2 + ys1 * k^2)

      # below
      xb <- c(xe2 * (1-k)^2 + xs2 * k^2)
      yb <- c(ye2 * (1-k)^2 + ys2 * k^2)

      # link
      xl <- c(xa, xs, rev(xb), rev(xe))
      yl <- c(ya, ys, rev(yb), rev(ye))

      xlink <- c(xlink, xl)
      ylink <- c(ylink, yl)
      set <- c(set, rep(i, length(xl)))
      group <- c(group, rep(df$trueGrp[i],
                            length(xl)))
    }
  }

  finalDat <- data.frame(xlink = xlink, ylink = ylink,
                         set = set, trueGrp = group,
                         stringsAsFactors = FALSE)
  return(finalDat)
}

# labels for the sets (or polygon items)
labelValue <- function(DF, r, ds , left = TRUE){

  mt <- c()
  trG <- c()
  angle <- c()
  x <- c()
  y <- c()
  sLab <- c()
  if(left){
    Gi <- DF$estGrp
    Si <- DF$isSig.end
    theta <- DF$theta.end
    labs <- DF$lab.end
  }else{
    Gi <- DF$trueGrp
    Si <- DF$isSig.start
    theta <- DF$theta.start
    labs <- DF$lab.start
  }

  for(i in unique(Gi)){

    sig.i <- unique(Si[Gi == i])
    lab.i <- unique(labs[Gi == i])

    if(sig.i){
      r.i <- r
    }else{r.i <- r - ds}

    ti <- theta[Gi == i]
    tmmi <- c(min(ti), max(ti))
    mti <- sum(tmmi)/2
    mt <- c(mt, mti)
    trG <- c(trG, i)
    y <- c(y, r.i * cos(mti))
    sLab <- c(sLab, lab.i)

    if(left){
      x <- c(x, -r.i * sin(mti))
      angle <- c(angle, mti/pi*180 + 270)
    }else{
      x <- c(x, r.i * sin(mti))
      angle <- c(angle, 90 - mti/pi*180)
    }

  }# i loop


  df <- data.frame(x = x, y = y, trueGrp = trG,
                   angle = angle, lab = sLab,
                   stringsAsFactors = FALSE)
  return(df)
}

# the outside circle : differential or non-differential
outSig <- function(df, rIn, rOut, itv, spPCT.r,
                   spPCT.l, gap.a, gap.b, ds){

  # Data for plot
  Num.l <- apply(df, 1, sum) # row -- left
  Num.r <- apply(df, 2, sum) # column -- right

  # left
  ppi <- pi - (gap.a + gap.b)

  if(length(Num.l)>1){
    sp.l <- ppi * spPCT.l / (length(Num.l)-1)
    }else{
    spPCT.l <- 0
    sp.l <- 0
  }
  space.l <- rep(sp.l, length(Num.l)-1)
  theta.l <- cumsum(Num.l)/sum(Num.l) * ppi * (1-spPCT.l)
  theta1.l <- c(0, head(theta.l, -1) + cumsum(space.l)) + gap.a
  theta2.l <- c(theta.l[1], tail(theta.l, -1) + cumsum(space.l)) +gap.a
  vSig.l <- names(Num.l) != "NoSignal"

  # ------------------ right -----------------------------
  if(length(Num.r)>1){
    sp.r <- ppi * spPCT.r / (length(Num.r)-1)
  }else{
    spPCT.r <- 0
    sp.r <- 0
  }
  space.r <- rep(sp.r, length(Num.r)-1)
  theta.r <- cumsum(Num.r)/sum(Num.r) * ppi * (1-spPCT.r)
  theta1.r <- c(0, head(theta.r, -1) + cumsum(space.r)) + gap.a
  theta2.r <- c(theta.r[1], tail(theta.r, -1) + cumsum(space.r)) + gap.a
  vSig.r <- names(Num.r) != "NoSignal"


  # if there are sets are inferred as significanltly differential
  if(sum(vSig.l) >0){
    ts.l <- c(itv[itv <= max(theta2.l[vSig.l]) & itv >= 0],
              max(theta2.l[vSig.l]))
    ts.r <- c(itv[itv <= max(theta2.r[vSig.r]) & itv >= 0],
              max(theta2.r[vSig.r]))
    x.s <- c(- rOut * sin(rev(ts.l)), rOut * sin(ts.r),
         rIn * sin(rev(ts.r)), - rIn * sin(ts.l),
         - rOut * sin(rev(ts.l))[1])

    y.s <- c(rOut * cos(rev(ts.l)), rOut * cos(ts.r),
         rIn * cos(rev(ts.r)), rIn * cos(ts.l),
         rOut * cos(rev(ts.l))[1])
  }else{ # if there are no significant sets in the inferred result
    ts.r <- c(min(theta1.r[vSig.r]),
              itv[itv <= max(theta2.r[vSig.r]) &
                    itv >= min(theta1.r[vSig.r])],
              max(theta2.r[vSig.r]))
    x.s <- c( rOut * sin(ts.r), rIn * sin(rev(ts.r)),
             rOut * sin(ts.r)[1])

    y.s <- c( rOut * cos(ts.r),
             rIn * cos(rev(ts.r)),
             rOut * cos(ts.r)[1])
  }

  # if there non-significant set in the inferred result
  if(sum(!vSig.l) >0){
    tn.l <- c(min(theta1.l[!vSig.l]),
              itv[itv >= min(theta1.l[!vSig.l]) & itv <= pi])
    tn.r <- c(min(theta1.r[!vSig.r]),
              itv[itv >= min(theta1.r[!vSig.r]) & itv <= pi])
    x.n <- c(- (rOut-ds) * sin(tn.l), (rOut-ds) * sin(rev(tn.r)),
           (rIn-ds) * sin(tn.r), - (rIn-ds) * sin(rev(tn.l)),
           - (rOut-ds) * sin(tn.l)[1])

    y.n <- c((rOut-ds) * cos(tn.l), (rOut-ds) * cos(rev(tn.r)),
           (rIn-ds) * cos(tn.r), (rIn-ds) * cos(rev(tn.l)),
           (rOut-ds) * cos(tn.l)[1])
  }else{ # if there no non-significant set in the inferred result
    tn.r <- c(min(theta1.r[!vSig.r]),
              itv[itv >= min(theta1.r[!vSig.r]) &
                    itv <= max(theta2.r[!vSig.r])],
              max(theta2.r[!vSig.r]))
    x.n <- c( (rOut-ds) * sin(rev(tn.r)),
             (rIn-ds) * sin(tn.r),
             (rOut-ds) * sin(rev(tn.r))[1])

    y.n <- c( (rOut-ds) * cos(rev(tn.r)),
             (rIn-ds) * cos(tn.r),
             (rOut-ds) * cos(rev(tn.r))[1])
  }

  df <- data.frame(x = c(x.s, x.n), y = c(y.s, y.n),
                   isSig = rep(c(TRUE, FALSE),c(length(x.s),
                                                length(x.n))),
                   stringsAsFactors = F)
  return(df)

  }

theme_blank <- ggplot2::theme(axis.line = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_blank(),
                              axis.text.y = ggplot2::element_blank(),
                              axis.ticks = ggplot2::element_blank(),
                              axis.title.x = ggplot2::element_blank(),
                              axis.title.y = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(),
                              panel.border = ggplot2::element_blank(),
                              panel.grid.major = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank(),
                              plot.background = ggplot2::element_blank())
