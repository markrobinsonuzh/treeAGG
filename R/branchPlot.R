#' Plot a tree branch
#'
#' \code{branchPlot} could show a tree plot of a selected branch with annotated
#' nodes and clades
#'
#' @param nodeBranch root node of a selected branch to be plotted (phylo object)
#' @param CladeLab a vector includes the labels of root node for interested clades
#'  and "other" (label for the branch part which are not included in the interested
#'  clades)
#' @param CladeCol a vector of colors to be assigned to the clades in CladeLab. The order between
#' CladeLab and CladeCol is matched
#' @param nodeLab a vector. Labels of nodes to be annotated in the plot; default is NULL,
#'  which indicates no node is annotated in the plot.
#' @param sizeNode a number to scale the size of the annotated node
#' @param add.Lab  indicate whether the label of annotated nodes should be added
#' @param LegendName the name of the legend title.
#'
#'
#' @return  a plot
#'
#' @example \notrun{
#' library(ape)
#' n <- 20
#' nodeLab <- paste("Node",1:(n-1),sep="")
#' tipLab <- paste("T",1:n,sep="")
#' # entire tree
#' Tree <- rtree(n,tip.label= tipLab)
#' Tree$node.label <- nodeLab
#'
#' # Subtrees
#' sTree <- subtrees(Tree)
#' names(sTree) <- nodeLab
#'
#' branchPlot(branch = sTree[["Node1"]],
#' CladeLab = c("Node17","Node12","other"),
#' CladeCol = c("blue","orange","black"),
#' nodeLab = c("Node14","Node13"),
#' LegendName = "Clade")
#' }
#'
#'
#'


branchPlot <- function(nodeBranch,stree,CladeLab,CladeCol,
                       nodeLab=NULL, sizeNode, add.Lab = TRUE,
                       LegendName, ...){

  p <- ggtree(stree[[nodeBranch]])
  names(CladeCol) <- CladeLab

  # reorder the CladeLab
  tl <- unlist(lapply(stree[CladeLab],
               FUN = function(x){length(x$tip.label)}))
  CladeLab1 <- CladeLab[order(tl,decreasing = TRUE)]
  CladeCol <- CladeCol[order(tl,decreasing = TRUE)]

  # group clades
  df <- p$data
  df[,"group"] <- 0
  df[,"colGroup"] <- tail(CladeCol,1)

  selClade <- head(CladeLab1, -1)
  # check whether selected clades exist in this branch
  ind.m <- match(selClade, df$label)
  if(sum(!is.na(ind.m)) ==0) {
    stop("None of selected clades are in this branch")
  }

  for (k in seq_len(length(selClade))){
    clade.k <- selClade[k]
    offsp <- FindOffspring(clade.k, stree,only.Tip = FALSE)
   idx <- match(offsp, df$label)
   df[idx,"group"] <- max(df[,"group"]) + 1
   df[idx,"colGroup"] <- CladeCol[clade.k]
   }

  df[,"group"] <- factor(df[,"group"])

  # indicated whether to annotate a node
   df[,"isLabel"] <- cladeData$label %in% nodeLab
   if(sum(df[,"isLabel"]) == 0){
    warning("None of the selected nodes are in this branch;
             Hence, no nodes would be annotated")
  }

   df1 <- df[,c("label","group","colGroup","isLabel")]

   new.CladeCol <- CladeCol[CladeCol %in% unique(df1$colGroup)]

   plot1 <- ggtree(stree[[branch]],aes(colour=colGroup),...) %<+% df1+
     scale_color_manual(name = LegendName,
                        values = unname(new.CladeCol),
                        labels = names(new.CladeCol))+
    theme(legend.position = "bottom")+
    guides(colour= guide_legend(override.aes = list(size = 2)),
           fill=FALSE, size=FALSE, label=FALSE)+
    geom_point2(aes(subset= isLabel,fill="blue"), size = sizeNode,
                color= "blue" ,show.legend = FALSE)

   if(add.Lab){
     plot1 +
       geom_text2(aes(subset= isLabel, label=label), color= "blue",
                  hjust=-0.3 ,show.legend = FALSE)
   }else{
     plot1
   }


}



