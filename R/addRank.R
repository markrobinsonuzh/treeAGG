#' create rank variable
#'
#' \code{addRank} is to create a rank variable which could be used later in
#' the tree aggregateion.
#'
#' @param ModOutp a dataframe output from differenetial analysis package (such as: edgeR)
#' @param varEst the name of the column which contains the main effect estimate
#' @param varSE  the name of the column which contains the standard error of main effect estimate
#' @param ...     see the arguments except data in \code{rvalues::rvalues}
#'
#' @return a dataframe
#' This function returns a new data frame which is actually the input dataframe ModOutp + a
#' new created column rvalue.
#'
#' rvalue is r value
#'
#' @examples {
#'
#' library(DESeq2)
#' # generate count tables from RNA-Seq data
#'  cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
#'  rownames(cnts) <- 1:nrow(cnts)
#'  cond <- factor(rep(1:2, each=5))
#'
#' # object construction
#' dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
#'
#' # standard analysis
#' dds <- DESeq(dds)
#' res <- results(dds)
#'
#' # compute r-value
#' Rres <- addRank(ModOutp = res, varEst = "log2FoldChange", varSE = "lfcSE",
#' varRank = "rvalue")
#' #' }


addRank <- function(ModOutp,varEst, varSE, ...){

  if(is.null(rownames(ModOutp))){
    stop("rownames of the count table ModOutp need to be specified")
  }

  # a data frame with the first column containing the main estimates
  # and the second column containing the nuisance terms (like SE)
  Est0 <- ModOutp[,c(varEst,varSE)]
  Est1 <- Est0[!is.na(Est0[,varEst]),]

  # use rvalues package to calculate the r values
  rRank0 <- rvalues::rvalues(Est1)$rvalues
  names(rRank0) <- rownames(Est1)

  # units with NA p-value have NA rvalues
  rRank <- rep(NA, nrow(Est0))
  names(rRank) <- rownames(Est0)
  rRank[names(rRank0)] <- rRank0

  # include rvalues as a new column
  df <- cbind.data.frame(ModOutp, rvalue = rRank)

  return(df)
}

