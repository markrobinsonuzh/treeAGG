#' create rank variable
#'
#' \code{addRank} is to create a rank variable which could be used later in
#' the tree aggregateion.
#'
#' @param ModOutp a dataframe output from differenetial analysis package (such as: edgeR)
#' @param varEst the name of the column which contains the main effect estimate
#' @param varSE  the name of the column which contains the standard error of main effect estimate
#' @param varRank the name to be used for the column of rank variable
#' @param pkg "edgeR" or "other"; if "edgeR", standard error of main effect estimates are calculated
#' within the function


addRank <- function(ModOutp,varEst, varSE, varRank,pkg = "edgeR"){
  if(pkg == "edgeR"){

  }else{
  # a data frame with the first column containing the main estimates
  # and the second column containing the nuisance terms (like SE)
  Est0 <- cbind.data.frame(Est = ModOutp[,varEst],
                          EstSE = ModOutp[,varSE])
  Est1 <- Est0[,!is.na(Est0[,varEst])]

  # use rvalues package to calculate the r values
  rRank0 <- rvalues::rvalues(Est1)$rvalues
  names(rRank0) <- rownames(Est1)

  # units with NA p-value have NA rvalues
  rRank <- rep(NA, nrow(Est0))
  names(rRank) <- rownames(Est0)
  rRank[names(rRank0)] <- rRank0
  }


  # include rvalues as a new column
  df <- cbind.data.frame(ModOutP, varRank = rRank)

  return(df)
}

