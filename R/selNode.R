#' Select a branch
#'
#' \code{selNode} selects a branch meeting the specified criteria in the number
#' of leaves and the count proportion.
#'
#' @param tree A phylo object
#' @param data A count table from the real data or a list output from
#'   \code{parEstimate}.
#' @param minTip the minimum number of leaves in the selected branch
#' @param maxTip The maximum number of leaves in the selected branch
#' @param minPr The minimum count proportion of the selected branch in a sample
#' @param maxPr The maximum count proportion of the selected branch in a sample
#' @param skip A character vector of node labels. These nodes are not the
#'   descendants or the ancestors of the selected branch.
#' @param all TRUE or FALSE. Default is FALSE. If FALSE, the branch node of a
#'   branch, which meet the requirements and has the minimum count proportion,
#'   is returned; otherwise branch nodes of all branches meeting the
#'   requirements are returned.
#'
#' @export
#'
#' @return The node whose descendant branch has the lowest proportion
#'
#' @examples{
#' data("cytofTree")
#' data("cytofCount")
#' set.seed(1)
#' library(ape)
#' tree <- as.phylo(cytofTree)
#' sN <- selNode(tree = tree, minTip = 20,
#' maxTip = 50, minPr = 0.02, maxPr = 0.03,
#' skip = NULL, data = cytofCount)
#' }
#'

selNode <- function(tree, data, minTip = 0, maxTip = Inf,
                    minPr = 0, maxPr = Inf,
                    skip = NULL, all = FALSE){

  ##------------ descendant tips ------------
  # descendant tips for each internal node

  # proportion of internal nodes
  leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
  nodI <- setdiff(tree$edge[, 1], leaf)
  nodI <- transNode(tree = tree, input = nodI)
  desI <- lapply(nodI, findOS, tree = tree,
                 only.Tip = FALSE, self.include = TRUE)
  tipI <- lapply(nodI, findOS, tree = tree,
                 only.Tip = TRUE, self.include = TRUE)
  names(tipI) <- nodI
  numI <- unlist(lapply(tipI, length))

  ##------------ node proportions -----------
  # tip proportions estimated from real data
  pars <- parEstimate(data = data)$pi
  names(pars) <- transNode(tree = tree, input = names(pars))

  # proportion for each node
  propList <- lapply(tipI, FUN = function(x){
    sum(pars[as.character(x)])
  })
  nodP <- unlist(propList)

  ##---------- sample ---------------
  tt <- cbind.data.frame(node = names(nodP),
                         proportion = nodP,
                         numTip = numI,
                         stringsAsFactors =FALSE)

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
    tipS <- lapply(skip, findOS, tree = tree,
                   only.Tip = TRUE, self.include = TRUE)
    tipS <- unlist(tipS)
    # take those without overlaps
    rmp <- sapply(st$node, FUN = function(x){
      tx <- findOS(ancestor = x, tree = tree, only.Tip = TRUE,
                   self.include = TRUE)
      ix <- intersect(tipS, tx)
      length(ix) == 0
    })

    new.st <- st[rmp, ]
  } else {
    new.st <- st
  }

  # return the one has the lowest proportion
  #ind <- which.min(abs(new.st$proportion - minPr))
  #final <- new.st[ind,]
  if (all) {
    final <- new.st
  } else {
    final <- new.st[which.min(new.st$proportion), ]
  }

  return(final)

  }
