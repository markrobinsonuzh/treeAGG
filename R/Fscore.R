
#' Calculate F1 score
#'
#' \code{fscore} calculates the F1 measure score
#'
#' @param truth a vector specifys the true class
#' @param estimate a vector specifys the estimated clusters
#'
#' @return a numeric value
#'
#' @examples {
#' set.seed(1)
#' truth <- sample(4:5, size = 20, replace = TRUE)
#' estimate <- sample(1:3, size = 20, replace = TRUE)
#'
#' # F score
#' fs <- fscore(truth, estimate)
#' }
#'



fscore <- function(truth, estimate){
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

  # F score
  # F1 = 2 * (precision*recall)/(precision+recall)
  prc <- index/csum
  rec <- index/rsum

  fscore <- 2 *(prc*rec)/(prc+rec)
  return(fscore)
}



