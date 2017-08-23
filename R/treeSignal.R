#' label true signal and found signal in a tree (the legend part is weird, fix this later)
#'
#' \code{treeSignal} returns a tree plot with true signal location and found signal
#' location labelled with different colors
#'
#' @param ModOutp  a dataframe output from differenetial analysis package (such as: edgeR)
#' @param TipDiff  differential abundant tips (charactor vector)
#' @param colSel   3 colors selected to label signal location (true, found, matched)
#' @param sizeSel  3 sizes selected to label signal location (true, found, matched)
#' @param pchSel   3 point types selected to label signal location (true, found, matched)
#' @param alpha    value for the color transparency of points
#' @param wtree    the entire phylogenetic tree which is based to do the aggregation
#' @param stree    a list of subtrees which are cut at the nodes of phylogenetic tree
#' @param P.lim    threshold for adjusted p-value (significant level)
#' @param VarCol   the column name of adjusted p-value
#' @param LegendName  the title of legend
#' @param ptree    the part of the tree to be plotted
#'
#' @return a tree
#'
#'@examples \notrun{
#' load("/Volumes/fiona/simu_500/R package/TreeInfer/data/Data_calcRateP.RData")
#' load("/Volumes/fiona/simu_500/R package/TreeInfer/data/tree.RData")
#'
#' edg1 <- edgeRes[[1]]
#'
#' treeSignal (ModOutp=edg1, TipDiff=c(TipC15,TipC19),
#' colSel=c("blue","firebrick","darkorchid")
#' ,pchSel = c(1,1,1),
#' alpha=1,wtree=Wtree, stree=Stree, sizeSel = c(2,1,2),
#' LegendName=" ",P.lim = 0.05,
#' VarCol = "FDR",ptree=wtree,size=0.1)
#'}
treeSignal <- function(ModOutp, TipDiff, colSel ,sizeSel, pchSel, alpha,
                       wtree, stree, P.lim = 0.05, VarCol = "FDR",
                       LegendName, ptree, ...){

  # Location of true signals
  trueSig <- signalFind(TipDiff,wtree,stree)

  # Location of signals found
  findSig <- find.p.min(wtree, ResTipNode= ModOutp, stree, P.lim, VarCol)

  # color for target points
  p <- ggtree::ggtree(ptree,branch.length = "none",...)

  names(pchSel) <- names(colSel) <- names(sizeSel) <-  c("True","Find","Match")
  pchPoint <- colPoint <- sizePoint <- rep(NA,nrow(p$data))
  names(pchPoint) <- names(colPoint) <- names(sizePoint) <- p$data$label

  pchPoint[intersect(setdiff(trueSig,findSig),p$data$label)] <- pchSel["True"]
  pchPoint[intersect(setdiff(findSig,trueSig),p$data$label)] <- pchSel["Find"]
  pchPoint[intersect(intersect(findSig,trueSig),p$data$label)] <- pchSel["Match"]


  colPoint[intersect(setdiff(trueSig,findSig),p$data$label)] <- colSel["True"]
  colPoint[intersect(setdiff(findSig,trueSig),p$data$label)] <- colSel["Find"]
  colPoint[intersect(intersect(findSig,trueSig),p$data$label)] <- colSel["Match"]

  sizePoint[intersect(setdiff(trueSig,findSig),p$data$label)] <- sizeSel["True"]
  sizePoint[intersect(setdiff(findSig,trueSig),p$data$label)] <- sizeSel["Find"]
  sizePoint[intersect(intersect(findSig,trueSig),p$data$label)] <- sizeSel["Match"]


  pointDef <- data.frame(label=p$data$label,
                             pchPoint = pchPoint,
                             colPoint = colPoint,
                             sizePoint = sizePoint,
                         stringsAsFactors = FALSE)
  p <- p %<+% pointDef

  # Tree plot

   p + ggtree::geom_point2(aes(subset= (p$data$label %in% c(trueSig,findSig)),
                    color = I(colPoint),size=sizePoint),
                alpha=alpha)+
    ggplot2::theme(legend.position="bottom",legend.key=element_blank())+
    ggplot2::guides(shape=FALSE,size=FALSE)+
    ggplot2::scale_color_manual(name=LegendName,
                       values=c("blue","firebrick","darkorchid"),
                       labels=names(colSel))
}

