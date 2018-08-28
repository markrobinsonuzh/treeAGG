#' Test for differential abundance: method 'treeAGG-DA-edgeR'
#'
#' Test differential abundance of entities using functions from the
#' \code{\link{edgeR}} (Robinson et al. 2010, \emph{Bioinformatics}; McCarthy et
#' al. 2012, \emph{Nucleic Acids Research}) to fit models and calculate
#' moderated test for each entity. We have used
#' \code{\link[edgeR]{estimateGLMRobustDisp}} to estimate the dispersion. The
#' statistical methods implemented in the \code{edgeR} package were originally
#' designed for the analysis of gene expression data such as RNA-sequencing
#' counts. Here, we apply these methods to counts that might be from microbes or
#' cells.
#'
#' The experimental design must be specified using a design matrix. The customized design matrix could be given by \code{design}.
#'
#' Normalization for samples is automatically performed by \code{edgeR} package.
#' More details about the calculation of normalization factor could be found
#' from \code{\link[edgeR]{calcNormFactors}}. A sample might include entities
#' corresponding to leaf nodes and internal nodes of tree. Only entities
#' corresponding to leaf nodes are used to calculate the library size of each
#' sample. The reason is that the abundance of an entity, corresponding to an
#' internal node, is calculated by taking sum of the abundance from its
#' descendant leaf nodes.
#'
#' @param obj A treeSummarizedExperiment object.
#' @param design A numeric matrix. It must be of full column rank. Defaults to
#' use all columns of \code{colData} to create design matrix. Note: Users should
#' check whether the default created design matrix is exactly what they want or
#' create their own design matrix using \code{\link[stats]{model.matrix}}.
#' @param contrast numeric vector or matrix specifying one or more contrasts of
#'   the linear model coefficients to be tested equal to zero. Number of rows
#'   must equal to the number of columns of design. If NULL, the last
#'   coefficient will be tested equal to zero.
#' @param normalize A logical value, TRUE or FALSE. The default is TRUE.
#' @param method Normalization method to be used. See
#'   \code{\link[edgeR]{calcNormFactors}} for more details.
#' @param prior.count average prior count to be added to observation to shrink
#'   the estimated log-fold-changes towards zero
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom edgeR DGEList calcNormFactors estimateGLMRobustDisp glmFit glmLRT
#'   topTags
#' @export
#' @return A treeSummarizedExperiment

runEdgeR <- function(obj, design = NULL, contrast = NULL,
                     normalize = TRUE, method = "TMM",
                     prior.count = 0.125){

    # calculate library size
    linkD <- linkData(obj)
    objL<- obj[linkD[linkD$isLeaf, "rowID"], ]
    tableL <- assays(objL)
    libSize <- lapply(tableL, FUN = function(x) {
        rs <- rowsum(x, group = rep(1, nrow(x)))
        return(rs[1, ])})

    # create DEGList
    tableA <- assays(obj)
    y <- lapply(seq_along(tableA), FUN = function(x) {
        yy <- DGEList(tableA[[x]], remove.zeros = TRUE)
        yy$samples$lib.size <- libSize[[x]]
        yy
    } )

    # normalisation
    if (normalize) {
        y <- lapply(seq_along(tableA), FUN = function(x) {
            calcNormFactors(y[[x]], method = method)})
    } else {
        y <- y
    }

    # default design matrix
    if (is.null(design)) {
        design <- designMatrix(data = colData(obj))
    }

    # estimate dispersion
    y <- lapply(y, estimateGLMRobustDisp, design = design)

    # build model
    fit <- lapply(y, glmFit, design = design, prior.count = prior.count)
    lrt <- lapply(fit, glmLRT, contrast = contrast)
    tt <- lapply(lrt, topTags, n = Inf, adjust.method = "BH",
                 sort.by = "none")
    final <- lapply(tt, FUN = function(x) { x$table })

    # output result to metadata
    outP <- lapply(seq_along(final), function(x) {
        # nodeLab and nodeNum columns from rowData
        #  dt <- rowData(obj)[, c("nodeLab", "nodeNum")]

        # find rows deleted
        idx <- as.numeric(rownames(final[[x]]))
        idc <- setdiff(seq_len(nrow(obj)), idx)

        # add and rearrange rows so that the output is also row-wise
        # corresponding to the assays data.
        df <- DataFrame(final[[x]])
        dm <- matrix(NA, nrow = length(idc), ncol = ncol(df),
                     dimnames = list(idc, colnames(df)))
        dfc <- DataFrame(dm)
        dfA <- rbind(df, dfc)
        dfA <- dfA[order(as.numeric(rownames(dfA))),]
        return(dfA)

    })

    outL <- treeSummarizedExperiment(assays = outP,
                                     rowData = rowData(obj),
                                     metadata = list(design = design,
                                                     contrast = contrast),
                                     tree = treeData(obj),
                                     linkData = linkD)
    return(outL)
}


#' Create design matrix
#'
#' \code{designMatrix} creates a design matrix by expanding factors to a set of
#' dummay variables and epanding interactions similarly.
#'
#'
#' \code{designMatrix} creates a design matrix using the data extracted from
#' \code{data}. \code{cols} specifies the columns to extract.
#' \code{designMatrix} is built on \code{\link[stats]{model.matrix}} with
#' \code{contrasts.arg = NULL}.
#'
#' @param data A \code{data.frame} or \code{DataFrame}.
#' @param cols A numeric vector. Specify columns to include in the design of
#' model matrix. Default is to include all columns.
#'
#' @importFrom stats model.matrix as.formula
#' @keywords internal
#' @return a matrix.

designMatrix <- function(data, cols = NULL) {

    stopifnot(class(data) %in% c("data.frame", "DataFrame"))

    # if cols is null, use all columns.
    if (is.null(cols)) {
        cols <- seq_len(ncol(data))
    }

    # create design matrix
    terms <- colnames(data)[cols]

    formula <- as.formula(paste("~", paste(terms, collapse = " + ")))

    design <- model.matrix(formula, data = data)

    return(design)
}

