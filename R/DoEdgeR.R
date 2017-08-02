


Redge<-function(nSam,countTab,isTip,isAnalyze,prior.count
                ,normalize=TRUE,outPut){
  # define conditions
  grp<-factor(rep(1:2, each=nSam))

  # correct sample size
  SampSize.c<-colSums(countTab[isTip,])
  y<-DGEList(countTab[isAnalyze,],group = grp,remove.zeros=TRUE)
  y$samples$lib.size<-SampSize.c

  # normalisation
  if(normalize){
    y <- calcNormFactors(y)
  }else{y<-y}

  # construct design matrix
  design<-model.matrix(~grp)
  # estimate dispersion
  y<-estimateDisp(y,design = design)
  fit<-glmFit(y,design = design)
  lrt<-glmLRT(fit)
  tt<-topTags(lrt,n = Inf,adjust.method="BH")

  # estimate the predictive log fold changes
  predlfc<-predFC(y,design,prior.count = prior.count)
  rownames(predlfc)<-rownames(y$counts)
  tt$table$predLFC<-predlfc[rownames(tt$table),2]

  # tagwise dispersion
  disp<-y$tagwise.dispersion
  names(disp)<-rownames(y$countTab)
  tt$table$tag.disp<-disp[rownames(tt$table)]
  colnames(tt$table)[colnames(tt$table)=="logFC"]<-"log2FC"
  colnames(tt$table)[colnames(tt$table)=="PValue"]<-"pvalue"

  if(outPut=="table"){return(tt$table)}
  if(outPut=="lrt"){return(lrt)}
  if(outPut=="y"){return(y)}
}
