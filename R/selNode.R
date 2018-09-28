#' Select branches
#'
#' \code{selNode} selects branches meeting the specified criteria in the number
#' of leaves and the count proportion.
#'
#' @param obj A leafSummarizedExperiment object.
#' @param tree A phylo object
#' @param data A count table from the real data or a list output from
#'   \code{parEstimate}.
#' @param minTip the minimum number of leaves in the selected branch
#' @param maxTip The maximum number of leaves in the selected branch
#' @param minPr The minimum count proportion of the selected branch in a sample.
#' A value between 0 and 1.
#' @param maxPr The maximum count proportion of the selected branch in a sample.
#' A value between 0 and 1.
#' @param skip A character vector of node labels. These nodes are not the
#'   descendants or the ancestors of the selected branch.
#' @param all TRUE or FALSE. Default is FALSE. If FALSE, the branch node of a
#'   branch, which meet the requirements and has the minimum count proportion,
#'   is returned; otherwise branch nodes of all branches meeting the
#'   requirements are returned.
#'
#' @export
#' @return The node whose descendant branch has the lowest proportion
#' @author Ruizhu Huang
#' @examples
#' set.seed(1)
#' data("tinyTree")
#' toyTable <- matrix(rnbinom(40, size = 1, mu = 10), nrow = 10)
#' colnames(toyTable) <- paste(rep(LETTERS[1:2], each = 2), rep(1:2, 2), sep = "_")
#' rownames(toyTable) <- tinyTree$tip.label
#'
#'
#' dat <- parEstimate(data = toyTable)
#' (out1 <- selNode(tree = tinyTree, data = dat, all = TRUE))
#' (out2 <- selNode(tree = tinyTree, data = dat, minTip = 4, maxTip = 9,
#' minPr = 0, maxPr = 0.8, all = TRUE))
#'
#' ## provide obj
#' lse1 <- leafSummarizedExperiment(tree = tinyTree,
#'                                 assays = list(toyTable))
#' out <- selNode(obj = lse1, all = TRUE)
#'

selNode <- function(obj = NULL, tree = NULL, data = NULL,
                    minTip = 0, maxTip = Inf,
                    minPr = 0, maxPr = 1,
                    skip = NULL, all = FALSE){

    # if the obj is provided with a leafSummarizedExperiment object
    # use it; otherwise use tree and data
    if (inherits(obj, "leafSummarizedExperiment")) {
        if (any(!is.null(tree), !is.null(data))) {
            stop("Set tree and data as NULL when obj is given. \n")
        }
        obj <- parEstimate(data = obj)
        tree <- metadata(obj)$tree
        data <- metadata(obj)$assays.par
    } else {
        tree <- tree
        data <- data
    }

    ##------------ descendant tips ------------
    # descendant tips for each internal node

    # proportion of internal nodes
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    nodI <- setdiff(tree$edge[, 1], leaf)
    nodLab <- transNode(tree = tree, input = nodI,
                        use.original = TRUE,
                        message = FALSE)
    nodLab_alias <- transNode(tree = tree, input = nodI,
                      use.original = FALSE,
                      message = FALSE)
    # descendant nodes
    desI <- lapply(nodI, findOS, tree = tree,
                   only.Tip = FALSE, self.include = TRUE)
    # descendant leaves
    tipI <- lapply(nodI, findOS, tree = tree,
                   only.Tip = TRUE, self.include = TRUE)
    names(tipI) <- nodLab_alias

    # number of descendant leaves
    numI <- unlist(lapply(tipI, length))

    ##------------ node proportions -----------
    # tip proportions estimated from real data
    pars <- parEstimate(data = data)$pi
    # a vector of node numbers with node labels as names
    vnum <- transNode(tree = tree, input = names(pars),
                      use.original = FALSE, message = FALSE)

    # proportion for each node
    propList <- lapply(tipI, FUN = function(x){
        sum(pars[names(vnum[vnum %in% x])])
    })
    nodP <- unlist(propList)

    ##---------- sample ---------------
    if (any(duplicated(nodLab))) {
        tt <- cbind.data.frame(nodeNum = transNode(tree = tree, input = names(nodP),
                                                   use.original = FALSE),
                               nodeLab = nodeLab,
                               nodeLab_alias = nodeLab_alias,
                               proportion = nodP,
                               numTip = numI,
                               stringsAsFactors =FALSE)
    } else {
        tt <- cbind.data.frame(nodeNum = transNode(tree = tree, input = names(nodP),
                                                   use.original = FALSE),
                               nodeLab = nodeLab,
                               proportion = nodP,
                               numTip = numI,
                               stringsAsFactors =FALSE)
    }


    if (maxPr < min(tt$proportion)) {
        stop("maxPr defined is even lower than the minimum value of
         node proportion", signif(min(tt$proportion),2), "\n")
    }
    # only consider nodes with enough tips and
    # desired proportion level
    st <- tt[tt$numTip >= minTip &
                 tt$numTip <= maxTip &
                 tt$proportion >= minPr &
                 tt$proportion <= maxPr,]
    if (nrow(st) == 0) {
        stop("No nodes fullfill the requirements;
         try other settings
         for tip numbers or proportions")
    }
    # remove those overlapped
    if (!is.null(skip)) {
        if (is(skip, character)) {
            skip <- transNode(tree = tree, input = skip, use.original = FALSE,
                              message = FALSE)
        }
        tipS <- lapply(skip, findOS, tree = tree,
                       only.Tip = TRUE, self.include = TRUE,
                       return = "number")
        tipS <- unlist(tipS)

        rmp <- vapply(st$nodeNum, FUN = function(x){
            tx <- findOS(ancestor = x, tree = tree, only.Tip = TRUE,
                         self.include = TRUE, return = "number")
            ix <- intersect(tipS, tx)
            length(ix) == 0
        }, FUN.VALUE = TRUE)

        new.st <- st[rmp, ]
    } else {
        new.st <- st
    }

    # return the one has the lowest proportion if all = FALSE
    #ind <- which.min(abs(new.st$proportion - minPr))
    #final <- new.st[ind,]
    if (all) {
        final <- new.st
    } else {
        final <- new.st[which.min(new.st$proportion), ]
    }

    return(final)

}
