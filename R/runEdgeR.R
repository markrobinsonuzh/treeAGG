#' Test for differential abundance: method 'treeAGG-DA-edgeR'
#'
#' Test differential abundance of entities using functions from the
#' \code{\link{edgeR}} (Robinson et al. 2010, \emph{Bioinformatics}; McCarthy et
#' al. 2012, \emph{Nucleic Acids Research}) to fit models and calculate
#' moderated test for each entity. We have used
#' \code{\link[edgeR]{estimateDisp}} to estimate the dispersion. The
#' statistical methods implemented in the \code{edgeR} package were originally
#' designed for the analysis of gene expression data such as RNA-sequencing
#' counts. Here, we apply these methods to counts that might be from microbes or
#' cells.
#'
#' The experimental design must be specified using a design matrix. The
#' customized design matrix could be given by \code{design}.
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
#' @param contrast numeric vector specifying one contrast of
#'   the linear model coefficients to be tested equal to zero. Its length
#'   must equal to the number of columns of design. If NULL, the last
#'   coefficient will be tested equal to zero.
#' @param normalize A logical value, TRUE or FALSE. The default is TRUE.
#' @param method Normalization method to be used. See
#'   \code{\link[edgeR]{calcNormFactors}} for more details.
#' @param prior.count average prior count to be added to observation to shrink
#'   the estimated log-fold-changes towards zero. See \code{prior.count} in
#'   \code{\link[edgeR]{glmFit}}
#' @param use.assays A numeric vector. It specifies which matrix-like elements
#'   in assays will be used to do analysis.
#' @param adjust.method A character string stating the method used to adjust
#'   p-values for multiple testing, passed on to \code{\link[stats]{p.adjust}}.
#'   It could be "bonferroni", "holm", "hochberg", "hommel", "BH", or "BY".
#'
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#'   topTags
#' @export
#' @return A treeSummarizedExperiment
#' \item{assays}{A list of tables}
#' \item{rowData}{It stores the information of rows in \code{assays}, and the
#' tables extracted from a \code{DGELRT} object that is generated by
#' \code{\link[edgeR:glmFit]{glmLRT}}. The later is stored as the internal part
#' of the \code{rowData}. More details or example could be found in the vignette
#' \code{Example of data analysis}}
#' \item{colData}{NULL}
#' \item{metadata}{
#'    \itemize{
#'    \item \code{use.assays} which elements in the \code{assays} have been
#'    used to run differential abundance analysis.
#'    \item \code{design} the design matrix as input.
#'    \item \code{contrast} the contrast vector as input.
#'    \item \code{output_glmFit} the output from \code{\link[edgeR]{glmFit}}. A
#'    object of \code{\link[edgeR]{DGEGLM-class}}
#'    }
#' }
#' @examples
#'
#' library(S4Vectors)
#' set.seed(1)
#' y <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(y) <- paste(rep(LETTERS[1:3], each = 10), rep(1:10,3), sep = "_")
#' rownames(y) <- tinyTree$tip.label
#'
#' rowInf <- DataFrame(nodeLab = rownames(y),
#'                     var1 = sample(letters[1:3], 10, replace = TRUE),
#'                     var2 = sample(c(TRUE, FALSE), 10, replace = TRUE))
#' colInf <- DataFrame(gg = factor(sample(1:3, 30, replace = TRUE)),
#'                     group = rep(LETTERS[1:3], each = 10))
#' toy_lse <- leafSummarizedExperiment(tree = tinyTree, rowData = rowInf,
#'                                     colData = colInf,
#'                                     assays = list(y, (2*y), 3*y))
#'
#' toy_tse <- nodeValue(data = toy_lse, fun = sum, tree = tinyTree,
#' message = TRUE)
#'
#' # build the model
#' contrastList <- list(contrast1 = c(0, 0, 0, -1, 1),
#'                      contrast2 = c(0, -1, 1, 0, 0))
#' mod <- runEdgeR(obj = toy_tse, contrast = contrastList)
#' # show results gained from the second element of the assasy
#' # sort by PValue
#' topNodes(mod, sort.by = "PValue", use.assays = 2)

runEdgeR <- function(obj, design = NULL, contrast = NULL,
                     normalize = TRUE, method = "TMM",
                     adjust.method = "BH",
                     prior.count = 0.125,
                     use.assays = NULL){

    # which elements from assays will be used for analysis.
    if (is.null(use.assays)) {
        use.assays <- seq_along(assays(obj))
    } else {
        use.assays <- use.assays
    }

    # calculate library size
    linkD <- linkData(obj)
    objL <- obj[linkD$isLeaf, ]
    tableL <- assays(objL, withDimnames = FALSE)[use.assays]
    libSize <- lapply(seq_along(tableL), FUN = function(x) {
        rs <- rowsum(tableL[[x]], group = rep(1, nrow(tableL[[x]])))
        rs[1, ]

    })

    # extract elements from assays for analysis
    tableA <- assays(obj, use.nodeLab = TRUE,
                     withDimnames = TRUE)[use.assays]

    # create DEGList
    y <- lapply(seq_along(tableA), FUN = function(x) {
        yy <- DGEList(tableA[[x]], remove.zeros = TRUE)
        yy$samples$lib.size <- libSize[[x]]
        yy
    } )

    # normalisation
    if (normalize) {
        y <- lapply(seq_along(tableA), FUN = function(x) {
            calcNormFactors(y[[x]], method = method)
        })
    } else {
        y <- y
    }

    # default design matrix
    if (is.null(design)) {
        design <- .designMatrix(data = colData(obj))
    } else { design <- design }

    # estimate dispersion
    y <- lapply(seq_along(y), FUN = function(x) {
        y <- estimateDisp(y[[x]], design = design)
        return(y)})

    # build model
    fit <- lapply(seq_along(y), FUN = function(x) {
        glmFit(y[[x]], design = design, prior.count = prior.count)})

    tt1 <- lapply(seq_along(fit), FUN = function(x) {
        # make sure a list is provided for the contrast
        if (is.null(contrast)) {
            lrt <- glmLRT(fit[[x]], contrast = contrast)
            tabx <- topTags(lrt, n = Inf, adjust.method = adjust.method,
                            sort.by = "none")$table
            xx <- list(tabx)
        } else {
            xx <- lapply(contrast, FUN = function(y) {
                lrt <- glmLRT(fit[[x]], contrast = y)
                tabx <- topTags(lrt, n = Inf,
                                adjust.method = adjust.method,
                                sort.by = "none")$table
                return(tabx)
            })
        }
        names(xx) <- names(contrast)
        return(xx) })

    tt2 <- lapply(seq_along(tt1), function(x) {
        # add and rearrange rows so that the output is also row-wise
        # corresponding to the assays data.
        res <- lapply(tt1[[x]], function(y) {
            # find rows deleted
            idx <- rownames(y)
            od <- linkD$nodeLab_alias
            if (is.null(od)) {
                od <- linkD$nodeLab
            }
            idc <- setdiff(od, idx)

            # add deleted rows in the table
            df <- DataFrame(y)
            dm <- matrix(NA, nrow = length(idc), ncol = ncol(df),
                         dimnames = list(idc, colnames(df)))
            dfc <- DataFrame(dm)
            dfA <- rbind(df, dfc)
            # sort rows
            dfA <- dfA[od,]
            return(dfA) })

        return(res)
        })

    # reshape: convert a list into a dataFrame
    tt3 <- lapply(seq_along(tt2), FUN = function(j) {
        x <- tt2[[j]]
        #dx <- x[[1]][, 0]
        dx <- x[[1]][, 0]
        for (i in seq_along(x)) {
            xi <- x[[i]]
            nxi <- names(x)[i]
            if (is.null(nxi)) { nxi <- "contrastNULL"}
            dx[[nxi]] <- xi
        }
        return(dx)
        })
    names(tt3) <- paste("result_assay", use.assays, sep = "" )

    # reshape to a dataFrame: one column for a result from one assay table
    tt4 <- tt3[[1]][, 0]
    for (i in seq_along(tt3)) {
        nxi <- names(tt3)[i]
        tt4[[nxi]] <- tt3[[i]]
    }
    rowData(obj)[["Results_internal_treeAGG"]] <- as(tt4, "internal_rowData")

    # contrast should have the same length as the sub-element of tt3
    if (is.null(contrast)) {
        contrast <- lapply(seq_along(tt3[[1]]), function(x) {NULL})
        names(contrast) <- "contrastNULL"
    }

    ## put the analysis result in the metadata
    metadata(obj)$output_glmFit <- fit
    metadata(obj)$use.assays <- use.assays
    metadata(obj)$design <- design
    metadata(obj)$contrast <- contrast


    # return(outL)
    return(obj)
}


#' Create design matrix
#'
#' \code{.designMatrix} creates a design matrix by expanding factors to a set of
#' dummay variables and epanding interactions similarly.
#'
#'
#' \code{.designMatrix} creates a design matrix using the data extracted from
#' \code{data}. \code{cols} specifies the columns to extract.
#' \code{.designMatrix} is built on \code{\link[stats]{model.matrix}} with
#' \code{contrasts.arg = NULL}.
#'
#' @param data A \code{data.frame} or \code{DataFrame}.
#' @param cols A numeric vector. Specify columns to include in the design of
#' model matrix. Default is to include all columns.
#'
#' @importFrom stats model.matrix as.formula
#' @keywords internal
#' @return a matrix.

.designMatrix <- function(data, cols = NULL) {

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

