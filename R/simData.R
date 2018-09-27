#'Simulate different scenarios of abundance change in entities
#'
#'\code{simData} simulates different abundance patterns for entities under
#'different conditions. These entities have their corresponding nodes on a tree.
#'More details about the simulated patterns could be found in the vignette via
#'\code{browseVignettes("treeAGG")}.
#'
#'@param tree A phylo object. Only use when \code{obj} is NULL.
#'@param data A matrix, representing a table of values, such as count, collected
#'  from real data. It has the entities corresponding to tree leaves in the row
#'  and samples in the column. Only use when \code{obj} is NULL.
#'@param obj A leafSummarizedExperiment object that includes a list of
#'  matrix-like elements, or a matrix-like element in assays, and a phylo object
#'  in metadata. In other words, \strong{obj} provides the same information
#'  given by \strong{tree} and \strong{data}.
#'@param scenario \dQuote{S1}, \dQuote{S2}, or \dQuote{S3} (see \bold{Details}).
#'  Default is \dQuote{S1}.
#'@param from.A,from.B The branch node labels of branches A and B for which the
#'  signal is swapped. Default, both are NULL. In simulation, we select two
#'  branches (A & B) to have differential abundance under different conditions.
#'  One could specify these two branches or let \code{doData} choose. (Note: If
#'  \code{from.A} is NULL, \code{from.B} is set to NULL).
#'@param minTip.A The minimum number of leaves in branch A
#'@param maxTip.A The maximum number of leaves in branch A
#'@param minTip.B The minimum number of leaves in branch B
#'@param maxTip.B The maximum number of leaves in branch B
#'@param minPr.A A numeric value selected from 0 to 1. The minimum abundance proportion
#'  of leaves in branch A
#'@param maxPr.A A numeric value selected from 0 to 1. The maximum abundance proportion
#'  of leaves in branch A
#'@param ratio A numeric value. The proportion ratio of branch B to branch A.
#'  This value is used to select branches(see \bold{Details}). If there are no
#'  branches having exactly this ratio, the pair with the value closest to
#'  \code{ratio} would be selected.
#'@param adjB a numeric value selected from 0 and 1 (only for \code{scenario} is
#'  \dQuote{S3}). Default is NULL. If NULL, branch A and the selected part of
#'  branch B swap their proportions. If a numeric value, e.g. 0.1, then the
#'  selected part of branch B decreases to its one tenth proportion and the
#'  decrease in branch B is added to branch A. For example, assume there are two
#'  experimental conditions (C1 & C2), branch A has 10 and branch B has 40 in
#'  C1. If adjB is set to 0.1, then in C2 branch B becomes 4 and branch A 46 so
#'  that the total proportion stays the same.
#'@param pct The percentage of leaves in branch B that have differential
#'  abundance under different conditions (only for scenario \dQuote{S3})
#'@param nSam A numeric vector of length 2, containing the sample size for two
#'  different conditions
#'@param mu,size The parameters of the Negative Binomial distribution. (see mu
#'  and size in \code{\link[stats]{rnbinom}}). Parameters used to generate the
#'  library size for each simulated sample.
#'@param n A numeric value to specify how many count tables would be generated
#'  with the same settings. Default is one and one count table would be obtained
#'  at the end. If above one, the output of \code{doData} is a list of matrices
#'  (count tables). This is useful, when one needs multiple simulations.
#'@param fun A function to derive the count at each internal node based on its
#'  descendant leaves, e.g. sum, mean. The argument of the function is a numeric
#'  vector with the counts of an internal node's descendant leaves.
#'@param seed a numeric value. Set seed to get reproducible results.
#'
#'@importFrom dirmult dirmult
#'@importFrom S4Vectors metadata
#'@importFrom SummarizedExperiment assays
#'@export
#'
#'@return a list of objects \item{FC}{the fold change of entities correspondint
#'  to the tree leaves.} \item{Count}{a list of count table or a count table.
#'  Entities on the row and samples in the column. Each count table includes
#'  entities corresponding to all nodes on the tree structure.}
#'  \item{Branch}{the information about two selected branches.} \describe{
#'  \item{A}{the branch node label of branch A} \item{B}{the branch node label
#'  of branch B} \item{ratio}{the count proportion ratio of branch B to branch
#'  A} \item{A_tips}{the number of leaves on branch A} \item{B_tips}{the number
#'  of leaves on branch B} \item{A_prop}{the count proportion of branch A (a
#'  value not above 1)} \item{B_prop}{the count proportion of branch B (the
#'  maximum is 1a value not above 1)} }
#'
#'@details \code{simData} simulates a count table for entities which are
#'  corresponding to the nodes of a tree. The entities are in rows and the
#'  samples from different groups or conditions are in columns. The library size
#'  of each sample is sampled from a Negative Binomial distribution with mean
#'  and size specified by the arguments \code{mu} and \code{size}. The counts of
#'  entities, which are located on the tree leaves, in the same sample are
#'  assumed to follow a Dirichlet-Multinomial distribution. The parameters for
#'  the Dirichlet-Multinomial distribution are estimated from a real data set
#'  specified by the argument \code{data} via the function \code{dirmult} (see
#'  \code{\link[dirmult]{dirmult}}). To generate different abundance patterns
#'  under different conditions, we provide three different scenarios,
#'  \dQuote{S1}, \dQuote{S2}, and \dQuote{S3} (specified via \code{scenario}).
#'  Our vignette provides figures to explain these three scenarios (try
#'  \code{browseVignettes("treeAGG")}). \itemize{ \item S1: two branches are
#'  selected to swap their proportions, and leaves on the same branch have the
#'  same fold change. \item S2: two branches are selected to swap their
#'  proportions. Leaves in the same branch have different fold changes but same
#'  direction (either increase or decrease). \item S3: two branches are
#'  selected. One branch has its proportion swapped with the proportion of some
#'  leaves from the other branch.}
#'
#'@author Ruizhu Huang
#'
#' @examples
#' \dontrun{
#' if(require(GUniFrac)){
#' data("throat.otu.tab")
#' data("throat.tree")
#'
#' # provide tree & data
#' count <- as.matrix(t(throat.otu.tab))
#' set.seed(1)
#' dat1 <- simData(tree = throat.tree, data = count, ratio = 2, seed = 1)
#'
#' # provide obj
#' treeDat <- treeSummarizedExperiment(tree = throat.tree,
#'                                     assays = list(count))
#' set.seed(1)
#' dat2 <- simData(obj = treeDat, ratio = 2, seed = 1123)
#' }
#'}


simData <- function(tree = NULL, data = NULL,
                    obj = NULL, scenario = "S1",
                    from.A = NULL, from.B = NULL,
                    minTip.A = 0, maxTip.A = Inf,
                    minTip.B = 0, maxTip.B = Inf,
                    minPr.A = 0, maxPr.A = 1,
                    ratio = 2, adjB = NULL,
                    pct = 0.6, nSam = c(50, 50),
                    mu = 10000, size = 50,
                    n = 1, fun = sum, seed = 1){
    set.seed(seed)
    # -------------------------------------------------------------------------
    # provide (tree & data)
    if(missing(obj)) {
        if(missing(tree) | missing(data)) {
            stop("tree or data is not provided")
        } else {
            obj <- doData(tree = tree, data = data, scenario = scenario,
                          from.A = from.A, from.B = from.B,
                          minTip.A = minTip.A, maxTip.A = maxTip.A,
                          minTip.B = minTip.B, maxTip.B = maxTip.B,
                          minPr.A = minPr.A, maxPr.A = maxPr.A,
                          ratio = ratio, adjB = adjB, pct = pct,
                          nSam = nSam, mu = mu, size = size,
                          n = n, fun = fun) }

    # -------------------------------------------------------------------------
    # provide obj
    } else {
        if(!inherits(obj, "leafSummarizedExperiment")){
            stop("obj should be a leafSummarizedExperiment object.")
        } else{
            # -------------------------------
            # don't use tree & data argument
            if ( (!missing(tree)) |
                 (!is.null(tree)) |
                (!missing(data)) |
                (!is.null(data)) ) {
                stop("Set tree = NULL and data = NULL when obj is a
                     leafSummarizedExperiment object. \n")
            }

            # confirme that the dirichlet multinomial parameters are available.
            # otherwise, estimate them.
            pars <- metadata(obj)$assays.par
            if (is.null(pars)) {
                obj <- parEstimate(data = obj)
                pars <- metadata(obj)$assays.par
            }
            # -------------------------------
            # data isn't provided, use obj assays data
            # if more than one table in assays, use the first one
            if(length(assays(obj)) > 1){
                message("\n more than one table provided in the assays;
                            only the first one would be used. \n")}
            data <- assays(obj)[[1]]
        }
        obj <- doData(tree = metadata(obj)$tree, data = pars,
                      scenario = scenario, from.A = from.A,
                      from.B = from.B,
                      minTip.A = minTip.A, maxTip.A = maxTip.A,
                      minTip.B = minTip.B, maxTip.B = maxTip.B,
                      minPr.A = minPr.A, maxPr.A = maxPr.A,
                      ratio = ratio, adjB = adjB, pct = pct,
                      nSam = nSam, mu = mu, size = size,
                      n = n, fun = fun)

    }
    return(obj)
}



#'Simulate a count table
#'
#'\code{doData} creates a count table for all nodes of a tree under two
#'different groups such that the tree would have different abundance patterns in
#'the different conditions.
#'
#'@param tree A phylo object
#'@param data A matrix, representing a count table from real data. It has the
#'  entities corresponding to tree leaves in the row and samples in the column.
#'@param scenario \dQuote{S1}, \dQuote{S2}, or \dQuote{S3} (see \bold{Details}).
#'  Default is \dQuote{S1}.
#'@param from.A,from.B The branch node labels of branches A and B for which the
#'  signal is swapped. Default, both are NULL. In simulation, we select two
#'  branches (A & B) to have differential abundance under different conditions.
#'  One could specify these two branches or let \code{doData} choose. (Note: If
#'  \code{from.A} is NULL, \code{from.B} is set to NULL).
#'@param minTip.A The minimum number of leaves in branch A
#'@param maxTip.A The maximum number of leaves in branch A
#'@param minTip.B The minimum number of leaves in branch B
#'@param maxTip.B The maximum number of leaves in branch B
#'@param minPr.A A numeric value selected from 0 to 1. The minimum abundance
#'  proportion of leaves in branch A
#'@param maxPr.A A numeric value selected from 0 to 1.The maximum abundance
#'  proportion of leaves in branch A
#'@param ratio The proportion ratio of branch B to branch A. This value is used
#'  to select branches(see \bold{Details}). If there are no branches having
#'  exactly this ratio, the pair with the value closest to \code{ratio} would be
#'  selected.
#'@param adjB a numeric value selected from 0 and 1 (only for \code{scenario} is
#'  \dQuote{S3}). Default is NULL. If NULL, branch A and branch B swap their
#'  proportions. If a numeric value, e.g. 0.1, then branch B decreases to its
#'  one tenth proportion and the decrease in branch B is added to branch A. For
#'  example, assume there are two experimental conditions (C1 & C2), branch A
#'  has 10 and branch B has 40 in C1. If adjB is set to 0.1, then in C2 branch B
#'  becomes 4 and branch A 46 so that the total proportion stays the same.
#'@param pct a numeric value selected from 0 and 1. The percentage of leaves in
#'  branch B that have differential abundance under different conditions (only
#'  for scenario \dQuote{S3})
#'@param nSam A numeric vector of length 2, containing the sample size for two
#'  different conditions
#'@param mu,size The parameters of the Negative Binomial distribution. (see mu
#'  and size in \code{\link[stats]{rnbinom}}). Parameters used to generate the
#'  library size for each simulated sample.
#'@param n A numeric value to specify how many count tables would be generated
#'  with the same settings. Default is one and one count table would be obtained
#'  at the end. If above one, the output of \code{doData} is a list of matrices
#'  (count tables). This is useful, when one needs multiple simulations.
#'@param fun A function to derive the count at each internal node based on its
#'  descendant leaves, e.g. sum, mean. The argument of the function is a numeric
#'  vector with the counts of an internal node's descendant leaves.
#'
#'@importFrom dirmult dirmult
#'@export
#'
#'@return a list of objects \item{FC}{the fold change of entities correspondint
#'  to the tree leaves.} \item{Count}{a list of count table or a count table.
#'  Entities on the row and samples in the column. Each count table includes
#'  entities corresponding to all nodes on the tree structure.}
#'  \item{Branch}{the information about two selected branches.} \describe{
#'  \item{A}{the branch node label of branch A} \item{B}{the branch node label
#'  of branch B} \item{ratio}{the count proportion ratio of branch B to branch
#'  A} \item{A_tips}{the number of leaves on branch A} \item{B_tips}{the number
#'  of leaves on branch B} \item{A_prop}{the count proportion of branch A (not
#'  above 1)} \item{B_prop}{the count proportion of branch B (not above 1)} }
#'
#'@details \code{doData} simulates a count table for entities which are
#'  corresponding to the nodes of a tree. The entities are in rows and the
#'  samples from different groups or conditions are in columns. The library size
#'  of each sample is sampled from a Negative Binomial distribution with mean
#'  and size specified by the arguments \code{mu} and \code{size}. The counts of
#'  entities, which are located on the tree leaves, in the same sample are
#'  assumed to follow a Dirichlet-Multinomial distribution. The parameters for
#'  the Dirichlet-Multinomial distribution are estimated from a real data set
#'  specified by the argument \code{data} via the function \code{dirmult} (see
#'  \code{\link[dirmult]{dirmult}}). To generate different abundance patterns
#'  under different conditions, we provide three different scenarios,
#'  \dQuote{S1}, \dQuote{S2}, and \dQuote{S3} (specified via \code{scenario}).
#'  \itemize{ \item S1: two branches are selected to swap their proportions, and
#'  leaves on the same branch have the same fold change. \item S2: two branches
#'  are selected to swap their proportions. Leaves in the same branch have
#'  different fold changes but same direction (either increase or decrease).
#'  \item S3: two branches are selected. One branch has its proportion swapped
#'  with the proportion of some leaves from the other branch.}
#'@author Ruizhu Huang
#'
#' @examples
#' \dontrun{
#' if(require(GUniFrac)){
#' data("throat.otu.tab")
#' data("throat.tree")
#'
#' dat <- doData(tree = throat.tree,
#' data = as.matrix(t(throat.otu.tab)),
#' ratio = 2)
#' }
#'}

doData <- function(tree, data, scenario = "S1",
                    from.A = NULL, from.B = NULL,
                    minTip.A = 0, maxTip.A = Inf,
                    minTip.B = 0, maxTip.B = Inf,
                    minPr.A = 0, maxPr.A = 1,
                    ratio = 2, adjB = NULL,
                    pct = 0.6, nSam = c(50, 50),
                    mu = 50, size = 10000,
                    n = 1, fun = sum){

    # ---check input is in correct format --------
    if(!inherits(tree, "phylo")){
        stop("tree should be a phylo object")
    }

    if (!inherits(data, "list")) {
        if (!inherits(data, "matrix")) {
            stop("data should be a matrix")
        } else {
            if (!setequal(rownames(data), tree$tip.label)) {
                stop("The rownames of data do not match with tree leaf labels")
            }
        }
    }


    # estimate parameters for Dirichlet-multinomial distribution
    data <- parEstimate(data = data)

    if (!is.null(from.A) && !is.null(from.B)) {
        pk <- infLoc(tree = tree, data = data,
        from.A = from.A, from.B = from.B)
    } else {
        pk <- pickLoc(tree = tree, data = data,
                      from.A  = from.A, minTip.A = minTip.A,
                      maxTip.A = maxTip.A, minTip.B = minTip.B,
                      maxTip.B = maxTip.B, minPr.A = minPr.A,
                      maxPr.A = maxPr.A, ratio = ratio)
    }


    beta <- doFC(tree = tree, data = data,
                 scenario = scenario,
                 branchA = pk$A,
                 branchB = pk$B,
                 ratio = pk$`ratio`,
                 adjB = adjB, pct = pct)



    count <- doCount(data = data, FC = beta,
                     nSam = nSam, mu, size, n)


    if(inherits(count, "list")) {
        grpDat <- data.frame(group = substr(colnames(count[[1]]), 1, 2))
        countLSE <- leafSummarizedExperiment(tree = tree, assays = count,
                                             metadata = list(
                                                 FC = beta,
                                                 branch = pk,
                                                 scenario = scenario),
                                             colData = grpDat)
    }

    if(inherits(count, "matrix")) {
        grpDat <- data.frame(group = substr(colnames(count), 1, 2))
        countLSE <- leafSummarizedExperiment(tree = tree,
                                             assays = list(count),
                                             metadata = list(
                                                 FC = beta,
                                                 branch = pk,
                                                 scenario = scenario),
                                             colData = grpDat)
    }
    obj <- nodeValue(data = countLSE, fun = sum)
    return(obj)
}


#' Select branches
#'
#' \code{pickLoc} selects two branches which meet the criteria specified by
#' the arguments
#'
#' @param tree A phylo object
#' @param data A count table (a matrix or a data frame). It has tree leaves in
#' rows and samples from different conditions in columns.
#' @param from.A The branch node label of branch A. In simulation, we select two
#' branches (A & B) to have differential abundance under different conditions.
#' If from.A is specified, then branch A is fixed. If from.A is NULL, one can
#' find a suitable branch which meets the criteria specified in \code{minTip.A},
#' \code{maxTip.A}, \code{minPr.A} and \code{maxPr.A}.
#' @param minTip.A The minimum number of leaves in branch A
#' @param maxTip.A The maximum number of leaves in branch A
#' @param minPr.A The minimum abundance proportion of leaves in branch A
#' @param maxPr.A The maximum abundance proportion of leaves in branch A
#' @param minTip.B The minimum number of leaves in branch B
#' @param maxTip.B The maximum number of leaves in branch B
#' @param ratio The proportion ratio of branch B to branch A. This value is
#' used to select branches(see \bold{Details}). If there are no branches having
#' exactly this ratio, the pair with the value closest to \code{ratio} would
#' be selected.
#'
#' @return a data frame of one row
#' @author Ruizhu Huang
#' @keywords internal

pickLoc <- function(tree, data, from.A,
                    minTip.A, maxTip.A,
                    minTip.B, maxTip.B,
                    minPr.A, maxPr.A, ratio) {

    # tip proportions estimated from real data
    pars <- parEstimate(data = data)$pi

    # proportion of internal nodes
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    nodI <- setdiff(tree$edge[, 1], leaf)
    desI <- lapply(nodI, findOS, tree = tree,
                   only.Tip = TRUE,
                   self.include = TRUE,
                   return = "label")
    nodP <- mapply(function(x, y) {
        sum(x[y])
    }, x = list(pars), y = desI)

    # matrix: abundance proprotion & the number of descendant leaves
    lenI <- unlist(lapply(desI, length))
    tt <- cbind(nodP, lenI)
    rownames(tt) <- transNode(tree = tree, input = nodI)

    # return error when the given limits for
    # proportion are outside
    # the range estimated from the real data.
    if (maxPr.A < min(tt[, 1])) {
        stop("maxPr.A is lower than the minimum value of
         node proportion", signif(min(tt[, 1]), 2), "\n")
    }
    if (minPr.A*ratio > max(tt[, 1])) {
        stop("minPr.A*ratio is above the maximum value of
         node proportion; try lower ratio", signif(min(tt[, 1]), 2), "\n")
    }

    # only consider nodes with enough tips and
    # desired proportion level
    if (is.null(from.A)) {
        tt.sel <- tt
    } else {
        # if node numbers, change them to node labels
        if (inherits(from.A, "character")) {
            from.A <- from.A
        } else {
            from.A <- transNode(tree = tree, input = from.A)
        }
        tt.sel <- tt[match(from.A, rownames(tt)), , drop = FALSE]
    }

    st <- tt.sel[tt.sel[, 2] >= minTip.A &
                     tt.sel[, 2] <= maxTip.A &
                     tt.sel[, 1] >= minPr.A &
                     tt.sel[, 1] <= maxPr.A, , drop = FALSE]
    if (nrow(st) == 0) {
        stop("No nodes fullfill the requirements;
         try other values for minTip.A, maxTip.A,
         minPr.A, or maxPr.A")
    }
    st2 <- tt[tt[, 2] >= minTip.B &
                  tt[, 2] <= maxTip.B, , drop = FALSE]

    # fold change between any two nodes (large/small)
    mm <- (1/st[, 1]) %o% st2[, 1]
    rownames(mm) <- rownames(st)

    maxM <- max(mm, na.rm = TRUE)
    minM <- min(mm, na.rm = TRUE)

    if (ratio >= 1 & ratio > maxM) {
        stop("could not find any two branches which fullfill
         these requirement;
         try lower ratio, lower minTip.A, or higher maxTip.B",
             "\n")
    }

    if (ratio <= 1 & ratio < minM) {
        stop("could not find any two branches which fullfill
         these requirement; try higher ratio or lower minTip.B",
             "\n")
    }

    # An entry in mm is set to NA if one node is the descendant of the other
    # nm <- matrix(sapply(seq_len(ncol(mm)), FUN = function(x) {
    #   # each column
    #   cn <- colnames(mm)
    #   cx <- cn[x]
    #
    #   # all rows
    #   rn <- rownames(mm)
    #   tx <- desI[rn]
    #
    #   cs <- lapply(tx,FUN = function(x) {
    #     length(intersect(x, desI[[cx]])) > 0
    #   })
    #   cv <- unlist(cs)
    #   fm <- mm[, x]
    #   fm[cv] <- NA
    #   fm}), nrow = nrow(mm))
    # preserve attributes of mm
    nm <- mm
    nm[] <- vapply(
        seq_len(ncol(mm)),
        FUN = function(x) {
            # each column
            cn <- colnames(mm)
            cx <- cn[x]

            # all rows
            rn <- rownames(mm)
            tx <- desI[rn]

            cs <- lapply(
                tx,
                FUN = function(x) {
                    length(intersect(x, desI[[cx]])) > 0
                }
            )
            cv <- unlist(cs)
            fm <- mm[, x]
            fm[cv] <- NA
            fm
        },
        FUN.VALUE = numeric(nrow(mm))
    )

    colnames(nm) <- colnames(mm)
    rownames(nm) <- rownames(mm)

    # iter <- 1
    # while (iter <= 200) {
    #   fs <- sample(rownames(nm),1)
    #   np <- nm[fs,]
    #
    #   ab <- abs(np-ratio)
    #   ind <- which( ab < 0.5 & ab==min(ab, na.rm = TRUE) )
    #   ind
    #   up <- np[ind[1]]
    #   up <- up[!is.na(up)]
    #   if( length(up) == 1 ){break}
    #   iter <-  iter+1
    # }
    # select the pair with the value closest to the ratio
    dif <- abs(nm - ratio)
    wi <- which(dif == min(dif), arr.ind = TRUE)
    si <- sample(seq_len(nrow(wi)), 1)
    an <- rownames(nm)[wi[si, 1]]
    bn <- colnames(nm)[wi[si, 2]]

    du <- cbind.data.frame(
        "A" = an,
        "B" = bn,
        "ratio" = round(nm[wi], digits = 2),
        "A_tips" = tt[an, 2],
        "B_tips" = tt[bn, 2],
        "A_prop" = round(tt[an, 1],
                            digits = 4),
        "B_prop" = round(tt[bn, 1],
                            digits = 4),
        stringsAsFactors =  FALSE
    )

    rownames(du) <- NULL
    return(du)

}

#' Provide the information of two branches
#'
#' \code{infLoc} is to give information of two branches about the count
#' proportion and the number of leaves
#'
#' @param tree A phylo object
#' @param data A count table (a matrix or a data frame). It has tree leaves in
#' rows and samples from different conditions in  columns.
#' @param from.A,from.B The branch node labels of Branch A, B.
#' @return A data frame of one row
#' @author Ruizhu Huang
#' @keywords internal

infLoc <- function(tree, data, from.A,
                   from.B) {

    # tip proportions estimated from real data
    pars <- parEstimate(data = data)$pi

    # proportion of internal nodes
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    nodI <- setdiff(tree$edge[, 1], leaf)
    desI <- lapply(nodI, findOS, tree = tree,
                   only.Tip = TRUE,
                   self.include = TRUE,
                   return = "label")
    nodP <- mapply(function(x, y) {
        sum(x[y])
    }, x = list(pars), y = desI)

    # matrix: abundance proprotion & the number of descendant leaves
    lenI <- unlist(lapply(desI, length))
    tt <- cbind(nodP, lenI)
    rownames(tt) <- transNode(tree = tree, input = nodI)

    # if both branches are given
    labA <- ifelse(is.character(from.A), from.A,
                   transNode(tree = tree, input = from.A))
    labB <- ifelse(is.character(from.B), from.B,
                   transNode(tree = tree, input = from.B))
    rAB <- tt[labB, 1]/tt[labA, 1]
    du <- cbind.data.frame(
        "A" = from.A,
        "B" = from.B,
        "ratio" = round(rAB, digits = 2),
        "A_tips" = tt[labA, 2],
        "B_tips" = tt[labB, 2],
        "A_prop" = round(tt[labA, 1],
                            digits = 4),
        "B_prop" = round(tt[labB, 1],
                            digits = 4),
        stringsAsFactors =  FALSE)

    rownames(du) <- NULL
    return(du)
}

#' Generate the fold change
#'
#' \code{doFC} generates fold changes for different scenarios
#'
#' @param tree A phylo object
#' @param data The real data (count table)
#' @param scenario Scenarios (\dQuote{S1}, \dQuote{S2}, \dQuote{S3})
#' @param branchA The branch node label of branch A.
#' @param branchB The branh node label of branch B.
#' @param ratio The proportion ratio between \code{branchB} and \code{branchA}
#' (B/A)
#' @param adjB A numeric value between 0 and 1 (only for \code{scenario}
#' is \dQuote{S3}). Default is NULL. If NULL, branch A and branch B swap their
#' proportions. If a numeric value, e.g. 0.1, then branch B decreases to its
#' one tenth proportion and the decrease in branch B is added to branch A.
#' For example, assume there are two experimental conditions (C1 & C2), branch
#' A has 10 and branch B has 40 in C1. If adjB is set to 0.1, then in C2 branch
#' B becomes 4 and branch A 46 so that the total proportion stays the same.
#' @param pct The percentage (in number) of the leaves in branchA that will
#' swap with branchB.

#'
#' @importFrom stats runif
#' @return numeric vector
#' @author Ruizhu Huang
#' @keywords internal

doFC <- function(tree, data, scenario,
                 branchA, branchB,
                 ratio, adjB, pct) {

    # beta: fold change for tips
    tips <- tree$tip.label
    beta <- rep(1, length(tips))
    names(beta) <- tips

    # tips on two branches
    tipA <- findOS(tree = tree, ancestor = branchA,
                   only.Tip = TRUE, self.include = TRUE,
                   return = "label")

    nod <- findOS(tree = tree, ancestor = branchA,
                  only.Tip = FALSE, self.include = TRUE,
                  return = "label")
    nodeA <- setdiff(nod, tipA)
    desA <- lapply(nodeA, findOS, tree = tree,
                   only.Tip = TRUE,
                   self.include = TRUE, return = "label")

    # swap proportion of two branches: tips in the same branch
    # have the same fold change
    if (scenario == "S1") {
        tipB <- findOS(tree = tree, ancestor = branchB,
                       only.Tip = TRUE, self.include = TRUE)

        beta[tipA] <- ratio
        beta[tipB] <- 1/ratio
    }

    # swap proportion of two branches: tips in the same branch
    # have different fold changes but same direction (either
    # increase or decrease)
    if (scenario == "S2") {
        tipB <- findOS(tree = tree, ancestor = branchB,
                       only.Tip = TRUE, self.include = TRUE)

        # proportion on two branches
        propA <- sum(data$pi[tipA])
        propB <- sum(data$pi[tipB])

        # swap proportions on two branches and randomly assign a fold change
        # value to the the leaves on a branch (log fold change in the same branch
        # should have the same sign)


            a1 <- runif(length(tipA))
            sa <- sum(a1 * propA)
            a2 <- (propB - propA)/sa
            a3 <- a1 * a2 + 1
            beta[tipA] <- a3

            b1 <- runif(length(tipB))
            sb <- sum(b1 * propB)
            b2 <- (propA - propB)/sb
            b3 <- b1 * b2 + 1
            beta[tipB] <- b3
    }

    # distribute signal randomly in one branch and evenly in
    # another branch
    # tip proportions estimated from real data
    pars <- parEstimate(data = data)$pi
    if (scenario == "S3") {
        iter <- 1
        while (iter <= 200) {
            # select only some tips
            selA <- sample(tipA, ceiling(length(tipA)*pct))
            # to make the selected tips disperse evenly in the branch
            subA <- lapply(desA, FUN = function(x) {
                ix <- intersect(x, selA)
                length(ix)/length(x)
            })
            subA <- unlist(subA)
            sumA <- sum(pars[selA])

            # the abundance proportion of the selected tips
            # are roughly equal to its number proportion in the branch
            # avoid (select all low or high abundance tips in the branch)

            spr <- sumA/sum(pars[tipA])
            ind.pr <- spr <= (pct + 0.05) & spr >= (pct - 0.05)

            if(all(subA <= 0.6) & ind.pr){break}
            iter <- iter+1
        }
        # ==================================================
        # cancel the random selection using selNode
        # ==================================================
        # if(is.null(branchB)){
        #   BranchL <- selNode(tree = tree,
        #                      minTip =length(selA),
        #                      maxTip = Inf,
        #                      minPr = max(ratio)*sumA,
        #                      maxPr = max(ratio)*sumA*1.5,
        #                      skip = branchA,
        #                      data = data)$node
        # }else{
        BranchL = branchB
        #}

        if (length(BranchL) == 0) {
            stop("No suitable branches.
           Try another branchA or another max of ratio... \n")
        }
        tipL <-  findOS(tree = tree, ancestor = BranchL,
                        only.Tip = TRUE, self.include = TRUE,
                        return = "label")
        sumL <- sum(pars[tipL])

        # ratio : 2, 5, 7, 8
        # if(length(ratio)==1){
        #   warning("For scenario S3, if multiple fold changes would be tried later
        #   to compare their results, it would be better to specify ratio as
        #           a vector so that the random tips selected would stay the same")
        # }

        # decide beta
        # beta <- sapply(ratio, FUN = function(x, y){
        #   y[selA] <- x
        #   xl <- (sumL-(x-1)*sumA)/sumL
        #   y[tipL] <- xl
        #   return(y)
        # }, y=beta)
        # colnames(beta) <- ratio
        if(is.null(adjB)){
            beta[selA] <- sumL/sumA
            beta[tipL] <- sumA/sumL
        }else{
            if(!is.numeric(adjB)){
                stop("adjB should be numeric")
            }
            beta[tipL] <- adjB
            beta[selA] <- (sumL*(1-adjB)+sumA)/sumA
            # beta[selA] <- fcA
            # beta[tipL] <- (sumL-(fcA-1)*sumA)/sumL
        }

    }

    return(beta)
}


#' generate a count table
#'
#' \code{doCount} generates a count table given some available information.
#' Here, the information includes the fold change between two conditions,
#' the parameters of Negative Binomial distritubion (which the sample library
#' size follows), a data table from real data (to estimate the proportion of
#' each entity in the sample), and the number of samples in two different
#' groups or conditions.
#'
#' @param data A matrix or data frame, representing the count table from real
#' data
#' @param FC A numeric vector, representing the fold changes
#' @param nSam A vector of length two containing the number of samples in two
#' different groups or conditions.
#' @param mu,size The parameters of the Negative Binomial distribution (see mu
#' and size in \code{\link[stats]{rnbinom}}.)
#' @param n A numeric value, representing the number of count tables generated.
#' If above one, the output is a list of matrices.
#'
#' @importFrom dirmult rdirichlet
#' @importFrom stats rmultinom rnbinom
#' @return a matrix or a list of matrices
#' @author Ruizhu Huang
#' @keywords internal
#'
doCount <- function(data, FC, nSam, mu,
                    size, n) {
    # parameters
    pars <- parEstimate(data = data)
    theta <- pars$theta
    gplus <- (1 - theta) / theta

    # tip proportion
    pr <- pars$pi
    p.c1 <- pr
    p.c2 <- pr * FC[names(p.c1)]

    # parameters for dirichlet distribution
    g.c1 <- p.c1 * gplus
    g.c2 <- p.c2 * gplus

    resList <- lapply(seq_len(n), FUN = function(j) {
        # condition 1
        n1 <- nSam[1]
        Mp.c1 <- matrix(0, nrow = n1, ncol = length(g.c1))
        rownames(Mp.c1) <- paste("C1_", seq_len(n1), sep = "")
        colnames(Mp.c1) <- names(p.c1)
        Mobs.c1 <- Mp.c1
        nSeq1 <- rnbinom(n = n1, mu = mu, size = size)
        for (i in seq_len(n1)) {
            Mp.c1[i, ] <- rdirichlet(1, g.c1)[1, ]
            Mobs.c1[i, ] <- rmultinom(1, nSeq1[i], prob = Mp.c1[i, ])[, 1]

        }

        # condition 2
        n2 <- nSam[2]
        Mp.c2 <- matrix(0, nrow = n2, ncol = length(g.c2))
        rownames(Mp.c2) <- paste("C2_", seq_len(n2), sep = "")
        colnames(Mp.c2) <- names(p.c2)
        Mobs.c2 <- Mp.c2
        nSeq2 <- rnbinom(n = n2, mu = mu, size = size)
        for (i in seq_len(n2)) {
            Mp.c2[i, ] <- rdirichlet(1, g.c2)[1, ]
            Mobs.c2[i, ] <- rmultinom(1, nSeq2[i], prob = Mp.c2[i, ])[, 1]
        }

        cb <- t(rbind(Mobs.c1, Mobs.c2))

        return(cb)}
    )

    if (n == 1) {
        count <-  do.call(rbind, resList)
    } else {
        count <- resList
    }

    return(count)
}



