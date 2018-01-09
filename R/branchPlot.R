#' Plot a tree branch
#'
#' \code{branchPlot} could show a tree plot of a selected branch with annotated
#' nodes and clades
#'
#' @param nodeBranch root node of a selected branch to be plotted (phylo object)
#' @param CladeLab a vector includes the labels of root node for interested clades
#' @param CladeCol a vector of colors to be assigned to the clades specified in CladeLab. The order between
#' CladeLab and CladeCol is matched. If NULL as default, then colors will be choosen automatically.
#' @param nodeLab a vector. Labels of nodes to be annotated in the plot; default is NULL,
#'  which indicates no node is annotated in the plot.
#' @param freqSel default NULL. Otherwise, it would be a numeric vector to specify the frequency 
#' of inferred nodes
#' @param sizeNode a number to scale the size of the annotated node
#' @param add.Lab  indicate whether the tip label of the inferred result should be added
#' @param colNode  colour specified for the nodes annotated (nodeLab)
#' @param LegendName the name of the legend title.
#' @param title plot title
#' @param scale numeric value
#'
#' @return  a plot
#'
#' @examples \notrun{
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
#' p1<- branchPlot(nodeBranch = "Node1",
#' stree = sTree,
#' CladeLab = c("Node17","Node12"),
#' CladeCol = c("blue","orange"),
#' nodeLab = c("Node14","Node13"),
#' sizeNode = 2,
#' colNode = "red",
#' LegendName = "Clade",zoomNode = c("Node14"),
#' scale=0.5, branch.length="none")
#' 
#' p2<- branchPlot(nodeBranch = "Node1",
#' stree = sTree,
#' CladeLab = c("Node17","Node12"),
#' CladeCol = NULL,
#' nodeLab = c("Node14","Node13"),
#' freqSel = c(3, 1),
#' sizeNode = 5,
#' colNode = "red",
#' LegendName = "Clade", branch.length="none")
#' 
#' multiplot(p1, p2, ncol=2)
#' }
#'
#'
#'


branchPlot <- function(nodeBranch,stree,CladeLab = NULL,CladeCol = NULL,
                       nodeLab=NULL, freqSel = NULL, sizeNode = 2, add.Lab = TRUE,
                       colNode = "darkviolet", LegendName, title = NULL,
                       scale = NULL, zoomNode = NULL, outData = FALSE, showTipLab = FALSE,...){
 
  if(length(nodeLab) == 0){ nodeLab <- NULL}
  if(length(CladeLab) == 0){ CladeLab <- NULL}
  if(is.factor(nodeLab)){stop("nodeLab is a factor; a character or numeric vector is required")}
  
  # assign colors for the specified clades 
  cl <- colors(distinct = TRUE)
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7",cl)
  #cbPalette <- c("grey19", "#E69F00", "#56B4E9", "#009E73", 
   #              "#F0E442", "#0072B2", "#D55E00", "#CC79A7",cl)
  if(is.null(CladeCol)){
    CladeCol <- tail(cbPalette, -1)[seq_along(CladeLab)]
  }else{CladeCol <- CladeCol}
  CladeLab <- c("other",CladeLab)
  CladeCol <- c(cbPalette[1], CladeCol)
  names(CladeCol) <- CladeLab

  # reorder the CladeLab
  tl <- unlist(lapply(stree[CladeLab],
               FUN = function(x){length(x$tip.label)}))
  CladeLab1 <- CladeLab[order(tl,decreasing = TRUE)]
  
  # group clades
  p <- ggtree::ggtree(stree[[nodeBranch]])
  df <- p$data
  df[,"group"] <- 0
  df[,"colGroup"] <- tail(CladeLab1,1)
  
  # check whether selected clades exist in the branch
  # with nodeBranch as root node
  selClade <- head(CladeLab1, -1)
  ind.m <- match(selClade, df$label)
  if(sum(!is.na(ind.m)) ==0 & !all(CladeLab=="other")) {
    warning("None of selected clades are in this branch")
  }

  # branch color
  for (k in seq_len(length(selClade))){
    clade.k <- selClade[k]
    offsp <- FindOffspring(clade.k, stree,only.Tip = FALSE)
    idx <- match(offsp, df$label)
    if(sum(!is.na(idx))>0){
      df[idx,"group"] <- k
      df[idx,"colGroup"] <- clade.k
    }
    
    }
   df[,"group"] <- factor(df[,"group"])

  # node annotated and size
   df[,"isLabel"] <- df$label %in% nodeLab
   df[,"Freq"] <- 0
   if(is.null(freqSel)){
     ind1 <- match(nodeLab, df$label)
     ind2 <- ind1[!is.na(ind1)]
     if(sum(ind2) >0){
       df[ind2,"Freq"] <- 1
     }
     
   }else{
     if(is.null(names(freqSel))){
       names(freqSel) <- nodeLab
       warning("The argument freqSel wasn't named and took the argument nodeLab as the name" )
     }
     
     df[match(nodeLab, df$label),"Freq"] <- freqSel[nodeLab]
   }
   
   
   if(sum(df[,"isLabel"]) == 0 & !is.null(nodeLab)){
    warning("None of the selected nodes are in this branch;
             Hence, no nodes would be annotated")
  }

   df1 <- df[,c("label","group","colGroup","isLabel", "Freq")]

   new.CladeCol <- CladeCol[names(CladeCol) %in% unique(df1$colGroup)]
   
   plot1 <- ggtree::ggtree(stree[[nodeBranch]],aes(colour=colGroup),...) %<+% df1+
     ggplot2::scale_color_manual(name = LegendName,
                        values = new.CladeCol)+
    ggplot2::theme(legend.position = "bottom", legend.text = element_text(size=25),
                   legend.key.size = unit(0.8,"cm"), legend.key.height = unit(0.4,"cm"))+
    ggplot2::guides(colour= guide_legend(override.aes = list(size = 2)),
           fill=FALSE, size=FALSE, label=FALSE)+
    ggtree::geom_point2(aes(subset= isLabel, size = Freq),fill= colNode,
                color= colNode, show.legend = FALSE) +
    ggplot2::scale_size_continuous(range= c(1, sizeNode))
    
   # ====================zoom ===========================
   if(!is.null(zoomNode) | !is.null(scale)){
     
     nodZ <- p$data$node[match(zoomNode,plot1$data$label)]
     if(is.null(zoomNode)|identical(zoomNode,character(0))){
       plot1 <- plot1
     }else{
       i <- 1
       repeat{
         plot1 <- plot1%>% scaleClade(nodZ[i], scale = scale)
         i <- i+1
         if(i > length(zoomNode)){
           break
         }
       }
     }
   }
   
   # ===============Tip label===========================
   if(showTipLab & (!is.null(zoomNode) | !identical(zoomNode,character(0))) ){
     zTip <- unlist(FindOffspring(ancestor = zoomNode, stree = stree, only.Tip = TRUE))
     #plot1 <- plot1 + geom_text2(aes(subset = (label %in% zTip), label=label),hjust=-0.3)
     plot1 <- plot1 + geom_text2(aes(label=label),hjust=-0.3, data=function(x){x[x$label %in% zTip,]})
   }else{
     plot1 <- plot1
   }
   
   if(add.Lab){
     sigTip <- unlist(FindOffspring(ancestor = nodeLab, stree = stree, only.Tip = TRUE))
     #plot1 <- plot1 + geom_text2(aes(subset = (label %in% zTip), label=label),hjust=-0.3)
     plot1 <- plot1 + geom_text2(aes(label=label, angle=angle),hjust=-0.3,
                                 data=function(x){x[x$label %in% sigTip,]},
                                 color = colNode)  + ggplot2::ggtitle(title)
   }else{
      plot1 <- plot1  + ggplot2::ggtitle(title)
    }
   # ===============================================
   #if(add.Lab){
    # plot1 +
     #  ggtree::geom_text2(aes(subset= isLabel, label=label), color= "blue",
      #            hjust=-0.3 ,show.legend = FALSE) + ggplot2::ggtitle(title)
   # }else{
    # plot1 + ggplot2::ggtitle(title)
   #}
   
   if(outData){
     return(plot1$data)
   }else{
     return(plot1)
   }
   
}



