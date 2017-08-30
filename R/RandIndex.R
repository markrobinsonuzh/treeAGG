
#' Calculate (adjusted) Rand Index
#'
#' \code{randIndex} calculates the (adjusted) Rand Index value
#'
#' @param truth a vector specifys the true class
#' @param estimate a vector specifys the estimated clusters
#' @param adjust a logical value; FALSE or TRUE; if TRUE, adjusted Rand Index is calculated;
#' otherwise, Rand Index is calculated
#'
#' @return a numeric value
#'
#' @examples {
#' set.seed(1)
#' truth <- sample(4:5, size = 20, replace = TRUE)
#' estimate <- sample(1:3, size = 20, replace = TRUE)
#'
#' # Rand Index
#' ri <- randIndex(truth, estimate, adjust = FALSE)
#' ari <- randIndex(truth, estimate, adjust = TRUE)
#' }
#'
randIndex <- function(truth, estimate, adjust = FALSE){
  rlen <- length(table(truth))
  ru <- unique(truth)
  clen <- length(table(estimate))
  cu <- unique(estimate)

  mat <- matrix(NA, nrow = rlen, ncol = clen)
  for(i in seq_len(rlen)){
    for(j in seq_len(clen)){
      mat[i, j] <- sum(truth == ru[i] & estimate == cu[j])
    }
  }

  index <- sum(choose(mat, 2))
  rsum <- sum(choose(rowSums(mat), 2))
  csum <- sum(choose(colSums(mat), 2))
  osum <- choose(length(truth), 2)
  Eindex <- rsum*csum/choose(length(truth), 2)
  Mindex <- (rsum + csum)/2

  ari <- (index - Eindex)/(Mindex - Eindex)
  ri <- (2*index + osum - rsum - csum)/osum

  if(adjust){
    return(ari)
  }else{
    return(ri)
  }

}



