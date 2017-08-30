glmnet <- function(countTab, nSam, isTip, isAnalyze,
                   normalize=TRUE, method = "TMM",
                   cv.fit = TRUE, ...){

  # define y response, site 1 as 0; site 2 as 1.
  y <- rep(c(0,1), nSam)

  # normalisation
  if(normalize){
    nFact <- calcNormFactors(countTab[isTip,], method = method)
    mat <- countTab %*% diag(nFact)
    mat <- t(mat)
  }else{mat <- t(countTab)}


  # lasso regularisation
  if(cv.fit){
    cv.mod <- glmnet::cv.glmnet(x = mat[, isAnalyze], y, family="binomial",
                                ...)
    return(cv.mod)
  }else{
    mod <- glmnet::glmnet(x = mat[, isAnalyze], y, family="binomial",
                          ...)
    return(mod)
  }
}

outLasso <- function(cvfit, fit){

  lambda <- cvfit$lambda.min
  coefs <- glmnet::coef.glmnet(fit, s = lambda)

  vSig <- rownames(coefs)[coefs[,1] != 0]
  vSig1 <- setdiff(vSig, "(Intercept)")
  return(vSig1)

}
