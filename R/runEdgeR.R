#' Differential analysis using edgeR package
#'
#' \code{runEdgeR} performs differential analysis for the data using the edgeR
#' package
#'
#' @param data  A numeric matrix corresponding to the count table of entities
#' (eg. genes, microbes or cell clusters) in samples under different conditions.
#' Each row represents an entity and each column represents a sample.
#' @param nSam A numeric vector. Each element represents the number of samples
#' in one condition. For example, if we have n1 samples in condition 1 and n2
#' samples in condition 2, then \strong{nSam} is (n1,n2).
#' @param isTip A logical vector. It has equal length to the number of rows in
#' \strong{data}. If TRUE, the corresponding row contributes to the library
#' size of a sample; otherwise, it doesn't. A provided count table might
#' include rows for leaves and internal nodes. Only the counts from leaf nodes
#' are really measured and should be included to do library size normalisation.
#' @param isAnalyze A logical vector with length equal to the row number of the
#' data. This is to indicate the row used for differential analysis. Default is
#' to use all rows.
#' @param prior.count A numeric value; the average prior count added to each
#' observation to compute estimated coefficients for a NB glm in such a way
#' that the log-fold-changes are shrunk towards zero
#' (see \code{\link[edgeR]{predFC}}.
#' @param normalize A logical value; indicating whether to correct for
#' compositionality when estimating size factors. If TRUE, TMM (weighted
#' trimmed mean of M-values) is applied; otherwise, raw library size is used.
#' @param method See \strong{method} in \code{\link[edgeR]{calcNormFactors}}
#'
#' @export
#' @importFrom edgeR DGEList calcNormFactors estimateGLMRobustDisp glmFit
#' glmLRT topTags predFC
#' @importFrom stats model.matrix p.adjust pnorm
#'
#' @return A data frame containing the elements
#'  logFC, the log-abundance ratio, i.e. fold change, for each tag in the
#'  two groups being compared logCPM, the log-average concentration/abundance
#'  for each tag in the two groups being compared
#'  PValue, exact p-value for differential expression using the NB model
#'
#'  FDR, the p-value adjusted for multiple testing as found using p.adjust
#'  using the method BH.
#'
#'  predLFC, the predictive log2 fold changes (SEE \code{\link[edgeR]{predFC}})
#'
#'  estimate, the estimated coefficient logFC = estimate / log(2)
#'
#'  tag.disp, the tagwise dispersion
#'
#'  waldAP, the adjusted p-value from Wald test ('BH' method; Benjamini &
#'  Hochberg (1995))
#' @author Ruizhu Huang
#' @examples
#'
#' data("cytofCount")
#' mod <- runEdgeR(data = cytofCount, nSam = c(5, 5),
#' isTip = rep(TRUE, nrow(cytofCount)), prior.count = 0,
#' normalize = TRUE)

runEdgeR <- function(data, nSam, isTip,
                     isAnalyze = NULL,
                     prior.count, normalize = TRUE,
                     method = "TMM") {
    # use all rows to do analysis
    if (is.null(isAnalyze)) {
        isAnalyze <- rep(TRUE, nrow(data))
    } else{
        isAnalyze <- isAnalyze
    }

    # define conditions
    grp <- factor(rep(seq_len(2), nSam))

    # correct sample size
    SampSize.c <- colSums(data[isTip, ])
    y <- edgeR::DGEList(data[isAnalyze, ], group = grp, remove.zeros = TRUE)
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
