tree.path<-function(TabTN,P.lim,wtree,stree
                    ,nClus,cluster,ClusDiff
                    ,collapse,show.minP,show.all
                    ,RangeAll,varCir,LegendName
                    ,title,alpha=0.9
                    ,size.minP,whole
                    ,pch.cir=c(1,19),...){

  # Create new variables and update (TabTN)
  TabTN.NEW<-NewVar.create(TabTN,P.lim,wtree,stree,show.minP)

  # All clusters
  for (i in 1:nClus){
    assign(paste("TipC",i,sep=""),names(which(cluster == i)))
  }

  # branch with differentially abundant tips
  # find node: 1. which is the farthest from root
  #            2. include the diff. abundant cluster
  tips<-lapply(stree,FUN = function(x){x$tip.label})
  if(whole){
    TipDiff<-unlist(tips)
  }else{
    TipDiff<-get(paste("TipC",ClusDiff,sep=""))
  }

  lt<-lapply(tips,FUN = function(x){all(TipDiff %in% x)})
  ln<-names(which(unlist(lt)==TRUE) )
  Nt<-unlist(lapply(stree[ln],FUN =function(x){x$Ntip}))
  tNode<-names(which(Nt==min(Nt)))
  BranchDiff<-stree[[tNode]]

  # branch tree plot
  p<-ggtree(BranchDiff,...)

  # model output result for this branch
  BranchTab<-TabTN.NEW[p$data$label,]

  # ========= decide branches to be omitted ===================
  if(whole){
    vnd2<-branch.omit(nClus,cluster,BranchDiff,TipDiff=NULL
                      ,stree)
  }else{
    vnd2<-branch.omit(nClus,cluster,BranchDiff,TipDiff,stree)
  }

  # ========= collapse or not =========================
  if(collapse){
    p<-branch.collapse(vnd2,p)
  }

  # =========== define edge line type =================
  # color the tips as grey if they are not in TipDiff
  group.lty<-rep(1,BranchDiff$Nnode+BranchDiff$Ntip)
  names(group.lty)<-p$data$label
  TipNDiff<-setdiff(BranchDiff$tip.label,TipDiff)
  #
  NodeNDiff<-lapply(tips,FUN = function(x){all(x %in% TipNDiff)})
  NND<-names(which(unlist(NodeNDiff)==TRUE))
  ind<-which(p$data$label %in% c(TipNDiff,NND,names(vnd2)))
  group.lty[ind]<-3
  group.lty<-factor(group.lty,levels = c(1,3))

  # ============ update p =================================
  BranchTab<-cbind.data.frame(label=p$data$label,BranchTab
                              ,groupLTY=group.lty)
  p <- p %<+% BranchTab
  p <- p+aes(linetype=groupLTY)


  # define color for different combinations of indicators above
  # isPosLFC, isNegLFC, isDp, isNaDp
  col.cirNeg<-rgb(0,0,205/255,alpha)
  col.cirPos<-rgb(238/255,118/255,33/255,alpha)
  col.cirNsNeg<-rgb(0,0,205/255,alpha/2)
  col.cirNsPos<-rgb(238/255,118/255,33/255,alpha/2)

  if(show.minP){
    cir.col<-c("Pos & nonSig"=col.cirNsPos
               ,"Neg & nonSig"=col.cirNsNeg
               ,"Pos & Sig & minP"=col.cirPos
               ,"Pos & Sig"=col.cirPos
               ,"Neg & Sig & minP"=col.cirNeg
               ,"Neg & Sig"=col.cirNeg)

    sel.pch<- c("Pos & nonSig"=1
                ,"Neg & nonSig"=1
                ,"Pos & Sig & minP"=2
                ,"Pos & Sig"=1
                ,"Neg & Sig & minP"=2
                ,"Neg & Sig"=1)
    cir.pch<-pch.cir[sel.pch]
    names(cir.pch)<-names(sel.pch)
  }else{
    cir.col<-c("Pos & nonSig"=col.cirNsPos
               ,"Neg & nonSig"=col.cirNsNeg
               ,"Pos & Sig"=col.cirPos
               ,"Neg & Sig"=col.cirNeg)
    sel.pch<- c("Pos & nonSig"=1
                ,"Neg & nonSig"=1
                ,"Pos & Sig"=1
                ,"Neg & Sig"=1)
    cir.pch<-pch.cir[sel.pch]
    names(cir.pch)<-names(sel.pch)
  }



  # =======determine whether to show non Significant nodes/tips=====
  #' display all nodes & tips
  if (show.all){
    if(varCir=="nLFDR.trunc"){
      if(show.minP){
        p<-p+
          geom_point2(aes(size=nLFDR.trunc,color=I(cir.lab)
                          ,shape=cir.lab))
      }else{
        p<-p+
          geom_point2(aes(size=nLFDR.trunc,color=I(cir.lab)))
      }

    }else{
      if(varCir=="log2FC"){
        if(show.minP){
          p<-p+
            geom_point2(aes(size=abs(log2FC),color=I(cir.lab)
                            ,shape=cir.lab))
        }else{
          p<-p+
            geom_point2(aes(size=abs(log2FC),color=I(cir.lab)))
        }
      }
    }

    #' show significant nodes & tips
  }else{
    if(varCir=="nLFDR.trunc"){
      if(show.minP){
        p<-p+ geom_point2(aes(size=nLFDR.trunc,color=I(cir.lab)
                              ,subset=(isDp &
                                         (isPosLFC | isNegLFC))
                              ,shape=cir.lab))
      }else{
        p<-p+ geom_point2(aes(size=nLFDR.trunc,color=I(cir.lab)
                              ,subset=(isDp & (isPosLFC | isNegLFC))))
      }
    }else{
      if(varCir=="log2FC"){
        if(show.minP){
          p<-p+ geom_point2(aes(size=abs(LFC),color=I(cir.lab)
                                ,subset=(isDp &
                                           (isPosLFC | isNegLFC))
                                ,shape=cir.lab))
        }else{
          p<-p+ geom_point2(aes(size=abs(LFC),color=I(cir.lab)
                                ,subset=(isDp & (isPosLFC | isNegLFC))))
        }
      }
    }
  }

  # Use RangeAll to make sure the point sizes are comparable among plots
  # otherwise the point sizes are only comparable in the same plot
  if(RangeAll){
    RangeVar<-range(TabTN.NEW[,varCir],na.rm = TRUE)
  }else{
    RangeVar<-range(BranchTab[,varCir],na.rm = TRUE)}

  p<-p+
    geom_tiplab(aes(subset=(isTip & isPosLFC & !isDp))
                ,color=col.cirNsPos,hjust = -0.5)+
    geom_tiplab(aes(subset=(isTip & isPosLFC & isDp ))
                ,color=col.cirPos,hjust = -0.5)+
    geom_tiplab(aes(subset=(isTip & isNegLFC & !isDp))
                ,color=col.cirNsNeg,hjust = -0.5)+
    geom_tiplab(aes(subset=(isTip & isNegLFC & isDp ))
                ,color=col.cirNeg,hjust = -0.5)+
    geom_tiplab(aes(subset=(isTip & isNaDp))
                ,color="grey",hjust = -0.5)+
    scale_size_continuous(limits=RangeVar)+
    ggtitle(label=title)+
    theme(legend.position="bottom"
          ,legend.key=element_blank())+
    guides(colour=guide_legend(override.aes=list(size=4,linetype=0),ncol=2)
           ,size=FALSE,alpha=FALSE,linetype=FALSE)+
    scale_color_manual(name=LegendName,values=cir.col)+
    scale_shape_manual(name=LegendName,values=cir.pch)


  # solve issue of too long tip label
  xRange<-ggplot_build(p)$layout$panel_ranges[[1]]$x.range
  p<-p+coord_cartesian(xlim = (xRange+c(0,xRange[2]/15)))

  # output
  print(p)

}
