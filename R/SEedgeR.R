

# library(edgeR)

SEedgeR <- function(y)
{
  X <- y$design
  coef <- y$coef
  phi <- y$dispersion
  mu <- y$fitted
  mu <- mu+1e-6
  nr <- nrow(mu)
  vb <- coef
  for (i in seq(nr))
  {
    W <- diag((mu[i,]^-1+phi[i])^-1)
    xtwx <- t(X) %*% W %*% X
    xtwxInv <- solve(xtwx)
    vb[i,] <- diag(xtwxInv %*% xtwx %*% xtwxInv)
  }
  se <- sqrt(vb)
  return(se)
}
