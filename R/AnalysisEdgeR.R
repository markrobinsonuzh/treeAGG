#' Differential analysis using edgeR package
#'
#' \code{Redge} is to do differential analysis for the data using edgeR package
#'
#' @param countTab  a numeric matrix indicates the count table of units (eg. genes or OTUs)
#'  in different samples under different conditions. Each row represents a unit and each
#'  column represents a sample
#' @param nSam a vector of length two to indicate the sample size for different conditions;
#'  e.g. (n1,n2) n1 is the sample size in condition 1 and n2 is the sample size in condition 2
#' @param isTip a logical vector with length equal to the row number of the countTab; This is
#'  to indicate the rows which contribute to the library size of a sample. For example, a count
#'  table might include both tip counts and node counts, but only the tip counts are gained
#'  from sequencing technology. The node counts are gained from aggregating the count of its descendant tips.
#'  In such situation, only the sum of tip counts are used to do library size normalisation
#' @param isAnalyze a logical vector with length equal to the row number of the countTab. This
#'  indicates the row which are included to do differential analysis.
#' @param prior.count a numeric value; the average prior count added to each observation to
#'  compute estimated coefficients for a NB glm in such a way that the log-fold-changes are
#'  shrunk towards zero (see \code{edgeR::predFC}).
#' @param normalize a logical value; indicate whether do normalization to scale library size
#'  if TRUE, TMM( weighted trimmed mean of M-values is applied); otherwise, raw library size
#'  is used.
#'
#' @return a data frame containing the elements
#'  logFC, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared
#'
#'  logCPM, the log-average concentration/abundance for each tag in the two groups being compared
#'
#'  PValue, exact p-value for differential expression using the NB model
#'
#'  FDR, the p-value adjusted for multiple testing as found using p.adjust using the method BH.
#'
#'  predLFC, the predictive log2 fold changes (SEE \code{edgeR::predFC})
#'
#'  estimate, the estimated coefficient logFC = estimate / log(2)
#'
#'  tag.disp, the tagwise dispersion
#'
#'  waldAP, the adjusted p-value from Wald test ('BH' method; Benjamini & Hochberg (1995))
#'  std.err the
#'
#' @examples {
#' library(dirmult)
#' library('GUniFrac')
#' data(throat.tree)
#' data(throat.otu.tab)
#'
#' # tree
#' Lab <- paste('Node',1:throat.tree$Nnode,sep='')
#' Wtree <- addNodeLab(treeO = throat.tree, nodeLab = Lab)
#' Stree <- pruneTree(wtree = Wtree)
#'
#' # count table for tips
#' DF.tip <- simuCount(RealDat = throat.otu.tab,
#' wtree = Wtree,nClus = 40,nSam = 5,
#' muNB = 10000,sizeNB = 5,
#' swapClus = c('cluster1','cluster19'),
#' diffClus = NULL,FC=NULL)
#'
#' # count table for tips and nodes
#' DF <- nodeCount(tipTable = DF.tip, wtree = Wtree, stree = Stree)
#'
#' # differential analysis
#' isNode <- substring(rownames(DF), 1, 4) == 'Node'
#' resDF <- Redge(countTab = DF, nSam = c(5,5), isTip = !isNode,
#' isAnalyze = rep(TRUE, nrow(DF)), prior.count =5, normalize=TRUE)
#'  }
#'


Redge <- function(countTab, nSam, isTip, isAnalyze, prior.count, normalize = TRUE, 
    method = "TMM") {
    # define conditions
    grp <- factor(rep(1:2, nSam))
    
    # correct sample size
    SampSize.c <- colSums(countTab[isTip, ])
    y <- edgeR::DGEList(countTab[isAnalyze, ], group = grp, remove.zeros = TRUE)
    y$samples$lib.size <- SampSize.c
    
    # normalisation
    if (normalize) {
        y <- edgeR::calcNormFactors(y, method)
    } else {
        y <- y
    }
    
    # construct design matrix
    design <- model.matrix(~grp)
    # estimate dispersion
    y <- edgeR::estimateGLMRobustDisp(y, design = design)
    fit <- edgeR::glmFit(y, design = design)
    lrt <- edgeR::glmLRT(fit)
    tt <- edgeR::topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "none")
    
    # estimate the predictive log fold changes
    predlfc <- edgeR::predFC(y, design, prior.count = prior.count)
    rownames(predlfc) <- rownames(y$counts)
    tt$table$predLFC <- predlfc[rownames(tt$table), 2]
    
    # wald test
    X <- fit$design
    coef <- fit$coefficients
    phi <- fit$dispersion
    mu <- fit$fitted.values
    mu <- mu + 1e-06
    nr <- nrow(mu)
    vb <- coef
    for (i in seq(nr)) {
        W <- diag((mu[i, ]^-1 + phi[i])^-1)
        xtwx <- t(X) %*% W %*% X
        xtwxInv <- solve(xtwx)
        vb[i, ] <- diag(xtwxInv %*% xtwx %*% xtwxInv)
    }
    waldStat <- coef/sqrt(vb)
    waldP <- (1 - apply(abs(waldStat), 2, pnorm)) * 2
    waldP.1 <- waldP[rownames(tt$table), 2]
    tt$table$waldAP <- p.adjust(waldP.1, method = "BH")
    tt$table$std.err <- sqrt(vb)[rownames(tt$table), 2]
    tt$table$estimate <- coef[rownames(tt$table), 2]
    
    # tagwise dispersion
    disp <- fit$dispersion
    names(disp) <- rownames(fit$coefficients)
    tt$table$tag.disp <- disp[rownames(tt$table)]
    # colnames(tt$table)[colnames(tt$table)=='logFC']<-'log2FC'
    # colnames(tt$table)[colnames(tt$table)=='PValue']<-'pvalue'
    
    return(tt$table)
}
