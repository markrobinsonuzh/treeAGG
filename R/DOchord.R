#' circular chord plot
#'
#' \{DOchord} creates a circular chord figure to compare the inferred result
#' and the truth. The truth is plotted in the right and the estimated (inferred)
#' is in the left. The two parts are separated by a middle black axis.
#'
#' @param infN a vector of location (tips/nodes) inferred to be
#'             significantly differential
#' @param trueN a vector of location (tips/nodes) which are really
#'              differentially abundant
#' @param stree a list of phylo object
#' @param decreasing a logical value with default FALSE. If FALSE,
#'                   the set (a set is a polygon in the figure) are
#'                   ordered increasingly so that the large sets are
#'                   arranged on the top (the size comparision is
#'                   made within sets w/wo signal).
#' @param spPCT.r  a numeric value between 0 and 1. The percentage of the
#'                 space in the right part of the figure. A high value
#'                 would lead to a wide gap between polygons in the figure
#' @param spPCT.l  the same to spPCT.r. The value is specified for the space
#'                 in left part of the figure
#' @param rOut a numeric value; the radius of the outside circle for polygon
#'             without signal (NoSignal)
#' @param rIn a numeric value; the radius of the inside circle for polygon
#'             without signal (NoSignal)
#' @param ds a numeric value; the difference between the radius of the circle
#'           for polygon with signal and that without signal
#' @param gap a numeric value between 0 and pi/2; the shortest distance between the
#'            middle axis and the polygon
#' @param size.labR a numeric value specified for the label size of polygons
#'                  in the right part of the figure
#' @param size.labL a numeric value specified for the label size of polygons
#'                  in the left part of the figure
#' @param alpha.link a numeric value between 0 and 1; specify the transparancy
#'                  of the link colours
#' @param color.linkb a color specified for the link boarder
#' @param size.linkb a numeric value specified for the link boarder size
#' @param color.group a vector specify the colours for the polygons in the
#'                    right part of the figure
#' @param color.sig a vector specify colours for Differential (abundance) and
#'                  Non-differential (abundance)
#' @param RI a numeric to be labelled on the plot to show the similarity
#'           between estimated result and the truth



DOchord <- function(infN, trueN, stree, decreasing = FALSE,
                    spPCT.r = 0.1, spPCT.l =0.1,
                    rOut = 1.1, rIn = 1,ds = 0.2,
                    gap.a = pi/90, gap.b = 0,
                    size.labL =2.5, size.labR = 3,
                    alpha.link = 0.7, color.linkb = "gray",
                    size.linkb = 0.1,
                    color.group = c("coral","cornflowerblue","gray48"),
                    color.sig = c("TRUE"="lightseagreen",
                                  "FALSE"="lightblue4"),
                    RI){

  # Data for plot
  df <- choreData(infN, trueN, stree, decreasing)

  # small intervals for polygon
  itv <- seq(from = 0, to = 180, by = 1) * (pi/180)

  # area : truth
  df.trueI <- axisValue(df = df, rIn = rIn, rOut = rOut,
                       itv = itv, spPCT = spPCT.r ,
                       gap.a = gap.a, gap.b = gap.b,
                       left = FALSE, ds = ds )

  # area : estimated
  df.estI <- axisValue(df = df, rIn = rIn, rOut = rOut,
                      itv = itv, spPCT = spPCT.l,
                      gap.a = gap.a, gap.b = gap.b,
                      left = TRUE, ds = ds)
  df.estI$x <- - df.estI$x

  # connect both areas
  # data : link end
  ke <- linkEnd(df = df, spPCT.l = spPCT.l, gap.a = gap.a, gap.b = gap.b)
  # data : link start
  ks <- linkStart(df = df, spPCT.r = spPCT.r, gap.a = gap.a, gap.b = gap.b)
  # data : end & start

  if(all(ke$trueGrp == ks$trueGrp)){
    df.both <- cbind.data.frame(ke[ ,colnames(ke)!= "trueGrp"], ks,
                                stringsAsFactors = FALSE)
  }else{stop("The link end and link start are not in the same order.")}

  df.link <- link(df = df.both, rR = rIn, ds = ds, rL = rOut)

  # significant or not
  df.out <- outSig(df= df, rIn=rOut+0.01, rOut = rOut+0.02,
                   itv, spPCT.r, spPCT.l, gap.a = gap.a,
                   gap.b = gap.b, ds)

  # add Labels
  df.labL <- labelValue(DF = df.both, r = rOut + 0.1, ds = ds, left = TRUE)
  df.labR <- labelValue(DF = df.both, r = rOut + 0.11, ds = ds, left = FALSE)

  hy.1 <- min(df.trueI$y[df.trueI$isSig])*((rOut+0.02)/(rIn))
  hx.1 <- sqrt((rOut+0.02)^2 - hy.1^2)
  hy.2 <- max(df.trueI$y[!df.trueI$isSig])*((rOut-ds+0.02)/(rOut-ds))
  hx.2 <- sqrt((rOut-ds+0.02)^2 - hy.2^2)


  df.labM <- data.frame(x = c(rIn/2, - rIn/2, hx.1+0.3, hx.2+0.3, rOut-1.5*ds),
                        y = c(0, 0, hy.1, hy.2, - rOut+1.5*ds ),
                        lab = c("Real", "Estimated", "Differential",
                                "Non-differential",
                                paste("RI : ", RI, sep = " ")),
                        colGrp = c("black", "black",
                                   color.sig[names(color.sig)=="TRUE"],
                                   color.sig[names(color.sig)=="FALSE"],
                                   "black"))

  # colors for polygon (right)
  names(color.group) <- unique(df.trueI$set)
  color.scale <- c(color.group, color.sig)

  ggplot2::ggplot()+
    ggplot2::geom_segment(aes(x = hx.1, y =  hy.1,
                     xend = hx.1+0.3, yend = hy.1),
                  color = color.sig[names(color.sig)=="TRUE"],
                  linetype = "dashed")+
    ggplot2::geom_segment(aes(x = hx.2, y =  hy.2,
                              xend = hx.2+0.4, yend = hy.2),
                          color = color.sig[names(color.sig)=="FALSE"],
                          linetype = "dashed") +
    ggplot2::geom_polygon(data = df.trueI,
                 aes(x, y, group = set, fill = factor(set)),
                 alpha = 0.7)+
    ggplot2::geom_polygon(data = df.link,
                 aes(x = xlink, y = ylink, group = set,
                     fill = factor(trueGrp)),
                 color=color.linkb,size=size.linkb,
                 alpha = alpha.link,
                 show.legend = F)+
    ggplot2::geom_polygon(data = df.estI,
                          aes(x, y, group = set),
                          fill= NA,color='black', size = 0.5)+
    ggplot2::geom_polygon(data = df.trueI, aes(x, y, group = set),
                          fill=NA,color='black', size = 0.5)+
    ggplot2::geom_polygon(data = df.out,
                          aes(x, y, fill = isSig),
                          show.legend = FALSE)+
    ggplot2::scale_fill_manual(values = color.scale)+
    ggplot2::geom_vline(xintercept = 0)+
    ggplot2::geom_label(data = df.labM, aes(x, y, label = lab,
                                  colour = I(colGrp)))+
    ggplot2::geom_text(data = df.labL,
              aes(x, y, label = lab,
                  angle = angle),
              size = size.labL)+
    ggplot2::geom_text(data = df.labR,
              aes(x, y, label = lab, angle = angle),
              size = size.labR)+
    ggplot2::theme(legend.position = "none")+
    theme_blank
}
