
# ======================
### class
# ======================

#  -----------------------------------------------------------------------------
#' internal_rowData: A virtual class
#' -----------------------------------------------------------------------------
#' A virtual class to be assigned to a column of \code{rowData} in a
#' \code{treeSummarizedExperiment} object. The columns with this class could be
#' optionally exported or hidden when using \code{rowData}.
#' @importFrom S4Vectors DataFrame
#' @keywords internal
setClass("internal_rowData", contains = "DataFrame")

#-------------------------------------------------------------------------------
#' phylo: A S3 class for the phylogenetic tree
#-------------------------------------------------------------------------------
#'
#' The \pkg{ape} package does not export its phylo class, probably because it
#' is not really defined formally anywhere. Technically, it is an S3 class
#' extended from the class list. Any exported definitions from the \pkg{ape}
#' package would be preferred to use if available.
#' @importFrom methods setOldClass
#' @keywords internal
phylo <- structure(list(), class = "phylo")
setOldClass("phylo")

#-------------------------------------------------------------------------------
# leafSummarizedExperiment
#-------------------------------------------------------------------------------
#' An S4 class leafSummarizedExperiment
#'
#' The class \strong{leafSummarizedExperiment} is an extension class of standard
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' class. It is designed to store rectangular data like a \code{matrix} for
#' entities (e.g., microbes or cell types), and the information of the
#' hiearchical structure of entities. The annotations on the rows and columns of
#' the rectangular data are strored in \code{rowData} and \code{colData},
#' respectively. Each row of the rectangular data could be mapped to a leaf node
#' of the tree via the rownames of \code{rowData} or a column named as
#' \code{nodeLab} in \code{rowData}. For example, the abundance count of
#' microbes collected from different samples could be stored in the
#' \code{assays} and the phylogenetic tree of the microbes in the
#' \code{metadata}. Each row of matrix-like data in \code{assays} represents one
#' microbial specie that could be mapped to a leaf node of the phylogenetic
#' tree. The link between the row and the node could be given via a column
#' (\code{nodeLab}), or the rownames of \code{rowData}. Each column of the
#' matrix-like data is a sample. The sample information is given in the
#' \code{colData}. The class \strong{leafSummarizedExperiment} has four slots
#' \code{assays}, \code{rowData} \code{colData} and \code{metadata} as the
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' class. More details about these four slots could be found in
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#'
#'
#' The \strong{leafSummarizedExperiment} class has more restrictions in data
#' structure than the \strong{SummarizedExperiment} class. In
#' \strong{leafSummarizedExperiment} class, it's required that
#' \itemize{
#' \item A \code{phylo} object is stored in \code{metadata} and named as
#' \code{tree}. The \code{phylo} object has a unique label for each leaf node.
#' \item Each rows of matrix-like data in \code{assays} could be mapped to a
#' leaf node of the tree. The label of leaf node is provided as the column
#' \code{nodeLab} in \code{rowData}, or as the rownames of \code{rowData} or
#' \code{assays}. If both provided, the information in column \code{nodeLab} is
#' used. }
#'
#' @section Constructor:
#' See \code{\link{leafSummarizedExperiment-constructor}} for constructor
#' functions.
#'
#' @section Accessor:
#' See \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} for
#' accessor functions.
#'
#' @importFrom methods setClass
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @name leafSummarizedExperiment-class
#' @exportClass leafSummarizedExperiment
#' @include classValid.R
#' @seealso \code{\link{treeSummarizedExperiment}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @author Ruizhu Huang
.lse <- setClass("leafSummarizedExperiment",
                 contains = "SummarizedExperiment",
                 validity = checkLSE)

#-------------------------------------------------------------------------------
# treeSummarizedExperiment
#-------------------------------------------------------------------------------
#' An S4 class treeSummarizedExperiment
#'
#' The class \strong{treeSummarizedExperiment} is an extension class of standard
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} class. It has six
#' slots. Four of them are traditional slots from
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} class:
#' \code{assays}, \code{rowData} \code{colData} and \code{metadata}. The other
#' two slots are \code{linkData} and \code{treeData}. The class
#' \strong{treeSummarizedExperiment} is designed to store rectangular data for
#' entities (e.g., microbes or cell types) (\code{assays}), information about
#' the hiearchical structure of entities (\code{treeData}), and information
#' about the mapping between the rows of the rectangular data and the nodes of
#' the tree (\code{linkData}).
#'
#' @slot linkData A \code{\link[S4Vectors]{DataFrame-class}} object. It gives
#'   map information between the rows of rectangular data and the nodes of tree.
#'   \itemize{
#'   \item \strong{nodeLab} The node labels on the tree.
#'   \item \strong{nodeLab_alias} An alias of column \code{nodeLab}. It is
#'   created only when there are missing value or duplicated value in column
#'   \code{nodeLab}. A prefix "Node_" and "Leaf_" is added to the node number
#'   (column \code{nodeNum}) for the internal nodes and the leaf nodes,
#'   respectively.
#'   \item \strong{nodeNum} The node numbers on the tree.
#'   \item \strong{isLeaf} This indicates whether a node is a leaf node.
#'   \item \strong{rowID} The row number in \code{assays}.
#'   }
#' @slot treeData A phylo object. It gives information about the hiearchical
#'   structure of the entities.
#' @slot ... See \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'   for more details about the slots inherited from \code{SummarizedExperiment}
#'   class.
#'
#' @section Constructor:
#' See \code{\link{treeSummarizedExperiment-constructor}} for constructor
#' functions.
#'
#' @section Accessor:
#' See \code{\link{treeSummarizedExperiment-accessor}} for accessor functions.
#'
#' @importFrom methods setClass
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @name treeSummarizedExperiment-class
#' @exportClass treeSummarizedExperiment
#' @seealso \code{\link{treeSummarizedExperiment}}
#'   \code{\link{treeSummarizedExperiment-accessor}}
#'   \code{\link{leafSummarizedExperiment}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
setClass("treeSummarizedExperiment",
         representation(linkData = "DataFrame",
                        treeData = "phylo"),
         contains = "SummarizedExperiment")




# ==========================================================================
### Constructor
# ==========================================================================
#' Construct a leafSummarizedExperiment object
#'
#' \code{leafSummarizedExperiment} is the constructor of the
#' \strong{leafSummarizedExperiment} class.
#'
#' @param tree A phylo object.
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#'
#' @importFrom utils head
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData
#'   "rowData<-" "metadata<-"
#' @importFrom S4Vectors metadata
#' @importFrom methods as
#' @export
#' @return A leafSummarizedExperiment object
#' @seealso \code{\link{leafSummarizedExperiment-class}}
#'   \code{\link{treeSummarizedExperiment-class}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @name leafSummarizedExperiment-constructor
#' @author Ruizhu HUANG
#' @examples
#' # tree
#' data("tinyTree")
#'
#' # the count table
#' count <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count) <- c(tinyTree$tip.label)
#' colnames(count) <- paste("C_", 1:10, sep = "_")
#'
#' # the sample information
#' sampC <- data.frame(condition = rep(c("control", "trt"),
#'                                     each = 5),
#'                     gender = sample(x = 1:2, size = 10,
#'                                     replace = TRUE))
#' rownames(sampC) <- colnames(count)
#'
#' # build a leafSummarizedExperiment object
#' lse <- leafSummarizedExperiment(tree = tinyTree,
#'                                 assays = list(count),
#'                                 colData = sampC)
#'
leafSummarizedExperiment <- function(tree, ...) {
    # -------------------------------------------------------------------------
    # create SummarizedExperiment object
    se <- SummarizedExperiment(...)

    # -------------------------------------------------------------------------
    # add tree in metadata
    metadata(se)$tree <- tree

    # -------------------------------------------------------------------------
    # keep only rows that could be assigned to tree leaves
    # use nodeLab if provided; otherwise, use rownames

    # use nodeLab
    tipLab <- tree$tip.label
    nodeLab<- rowData(se)$nodeLab
    if (!is.null(nodeLab)) {
        isIn <- (nodeLab %in% tipLab)
        isOut <- !isIn
        if (sum(isOut) > 0) {
            message(sum(isOut), "rows are removed from tables of *assays*. \n",
                nodeLab[isOut],
                " can't be matched to any leaf nodes of the tree. \n")}
        se <- se[isIn, ] }
    # use rownames
    if (is.null(nodeLab)) {
        rowNam <- rownames(se)
        if (is.null(rowNam)) {
            stop("Either the rownames of rowData or a nodeLab column in ",
                 "rowData should be provided \n.")
        }
        isIn <- rowNam %in% tipLab
        isOut <- !isIn
        if (sum(isOut) > 0) {
            message(sum(isOut), "rows are removed from tables of *assays*. \n",
                rowNam[isOut],
                " can't be matched to any leaf node of the tree. \n")}
        se <- se[isIn, ]
        }

    # output as leafSummarizedExperiment
    #lse <- as(se, "leafSummarizedExperiment")
    lse <- .lse(se)
    return(lse)
    }

#' Construct a treeSummarizedExperiment object
#'
#' \code{treeSummarizedExperiment} constructs a treeSummarizedExperiment object.
#'
#' @param tree A phylo object
#' @param linkData A data frame or \code{\link[S4Vectors]{DataFrame-class}}. It
#'   has the same number of rows as the matrix-like elements of \code{assays}.
#'   The row order of the \code{linkData} matches with that of  the matrix-like
#'   element of \code{assays}. It has the following columns.
#'   \itemize{
#'   \item \strong{nodeLab} The labels of nodes on the tree.
#'   \item \strong{nodeLab_alias} An alias of column \code{nodeLab}. It is
#'   created only when there are missing value or duplicated value in column
#'   \code{nodeLab}. A prefix "Node_" and "Leaf_" is added to the node number
#'   (column \code{nodeNum}) for the internal nodes and the leaf nodes,
#'   respectively.
#'   \item \strong{nodeNum} The numbers of nodes
#'   \item \strong{isLeaf} It indicates whether the node is a leaf node
#'   \item \strong{rowID} The corresponding row number of the matrix-like
#'   elements in \code{assays} }
#' @param message A logical value. Defaut is FALSE. If TRUE, the running process
#'   is given.
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
#' @export
#' @return a treeSummarizedExperiment object
#' @name treeSummarizedExperiment-constructor
#' @author Ruizhu HUANG
#' @seealso \code{\link{treeSummarizedExperiment-class}}
#'   \code{\link{treeSummarizedExperiment-accessor}}
#'   \code{\link{leafSummarizedExperiment-class}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @examples
#' data("tinyTree")
#'
#' # the count table
#' count <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count) <- c(tinyTree$tip.label)
#' colnames(count) <- paste("C_", 1:10, sep = "_")
#'
#' # The sample information
#' sampC <- data.frame(condition = rep(c("control", "trt"), each = 5),
#'                     gender = sample(x = 1:2, size = 10, replace = TRUE))
#' rownames(sampC) <- colnames(count)
#'
#' # build a treeSummarizedExperiment object
#' tse <- treeSummarizedExperiment(tree = tinyTree,
#'                                 assays = list(count),
#'                                 colData = sampC)
#'
#' # if there are data to be put in the linkData
#' tax <- data.frame(family = rep(LETTERS[1:2], each = 5))
#' tax$nodeLab <- tinyTree$tip.label
#' ntse <- treeSummarizedExperiment(tree = tinyTree,
#'                                 assays = list(count),
#'                                 colData = sampC,
#'                                 linkData = tax)
#'
treeSummarizedExperiment <- function(tree = NULL, linkData = NULL,
                                     message = FALSE, update.all = TRUE,
                                     ...){

    # -------------------------------------------------------------------------
    ## create a SummarizedExperiment object
    se <- SummarizedExperiment(...)

    # -------------------------------------------------------------------------
    ## tree: a phylo object
    if(!inherits(tree, "phylo")) {
        stop(tree, ": A phylo object is required", "\n")
    }

    # -------------------------------------------------------------------------
    ## create linkData: the number of rows equals to the number of nodes
    # the nodes
    if (message) {cat("Collect the information of nodes from the tree... \n")}

    emat <- tree$edge
    leaf <- setdiff(emat[, 2], emat[, 1])
    nodeI <- setdiff(emat[, 1], leaf)
    nodeA <- c(leaf, nodeI)

    ## Create the main part of the link data
    labA <- transNode(tree = tree, input = nodeA,
                      use.alias = FALSE, message = FALSE)
    if (anyDuplicated(labA)) {
        # provide the alias of node labels if there are duplicated values
        labA_a <- transNode(tree = tree, input = nodeA,
                            use.alias = TRUE, message = FALSE)
        linkD <- DataFrame(nodeNum = nodeA, nodeLab = labA,
                           nodeLab_alias = labA_a,
                           isLeaf = nodeA %in% leaf)
    } else {
        linkD <- DataFrame(nodeNum = nodeA, nodeLab = labA,
                           isLeaf = nodeA %in% leaf)
    }

    ## Create the extra part of the link data from the provided linkData
    if (!is.null(linkData)) {
        # message
        if (message) {cat("Combine the linkData provided... \n")}

        linkData <- DataFrame(lapply(linkData, as.character))
        labLK <- linkData$nodeLab_alias
        if (is.null(labLK)) {
            labLK <- linkData$nodeLab
        }
        if (is.null(labLK)) {
            stop("A column (nodeLab) or row names should be provided in linkData. \n ")
        }
        rownames(linkData) <- labLK

        # check wehther the linkData for the same node is the same
        numLK <- transNode(tree = tree, input = labLK,
                           use.alias = FALSE, message = FALSE)
        ind1 <- duplicated(linkData)
        ind2 <- duplicated(numLK)
        ind3 <- ind1 != ind2
        if (any(ind3)) {
            lab <- transNode(tree = tree, input = numLK[ind3],
                             use.alias = FALSE, message = FALSE)
            stop("The linkData provided for the nodes ",
                 lab, " are different. \n")
        }

        # remove duplicated data
        linkN <- linkData[, !colnames(linkData) %in% colnames(linkD),
                          drop = FALSE]
        linkN$numLK <- numLK
        linkN <- linkN[!duplicated(linkN), ]

        # include the given linkData as the extra part of link data
        linkC <- linkN[rep(1, nrow(linkD)), ]
        linkC[] <- NA
        ind4 <- match(linkD$nodeNum, linkN$numLK)
        linkC <- linkN[ind4, !names(linkN) %in% "numLK"]

        # update nodes based on the given linkData
        if (message) {cat("Update nodes that are not in the linkData... \n")}
        # if (update.all) {
        ind5 <- which(is.na(ind4))
        if (length(ind5)) {
            numOt <- linkD$nodeNum[ind5]
            funU <- function(x) {
                ux <- unique(x)
                ifelse(length(ux) == 1, ux, NA)
            }
            linkAD <- nodeValue.A(data = linkData, fun = funU,
                                  tree = tree, message = message,
                                  level = numOt)
            linkC[ind5, ] <- linkAD[, !names(linkAD) %in% "nodeLab"]
        }

        #}
        linkA <- cbind(linkD, linkC)

        } else {
        linkA <- linkD
    }

    # add rowID to store information about the map information between rows of
    # assays and nodes of tree
    # (1) if the nodeLab column exist, they should match with the labels of
    # tree nodes.
    if (message) {cat("Generate the rowID column in the link data... \n")}
    labRD <- rowData(se)$nodeLab
    if (!is.null(labRD)) {
        numRD <- transNode(tree = tree, input = labRD,
                           use.alias = FALSE, message = FALSE)
        # keep only rows that could be mapped to nodes of tree
        isIn <- numRD %in% linkA$nodeNum
        isOut <- !isIn
        if (sum(isOut) > 0) {
            message(sum(isOut), "rows are removed. They can't be matched to
                    any node of the tree. \n")}
        se <- se[isIn, ]
        rowID <- lapply(seq_along(linkA$nodeNum), FUN = function(x) {
            numRDL <- numRD[isIn]
            which(numRDL %in% x)

        })
        linkA$rowID <- rowID
        }

    # (2) if the nodeLab column doesn't exist, rownames should match with
    # the labels of tree leaves.
    if (is.null(labRD)) {
        # the row names
        labRN <- rownames(se)
        # if neither nodeLab nor rownames are provided, return error.
        if (is.null(labRN)) {
            stop("Either a nodeLab column or row names should be
                     provided for row data \n.")}
        # label to number
        numRN <- transNode(tree = tree, input = labRN,
                           use.alias = FALSE, message = FALSE)

        # keep only those could be mapped to nodes of the tree
        isIn <- numRN %in% linkA$nodeNum
        isOut <- !isIn
        if (sum(isOut) > 0) {
            message(sum(isOut), " rows of assays are removed.
                    They cannot be maped to any node of the tree. \n")}
        se <- se[isIn, ]
        rowID <- lapply(seq_along(linkA$nodeNum), FUN = function(x) {
            numRNL <- numRN[isIn]
            which(numRNL %in% x)
        })
        linkA$rowID <- rowID
    }

    # -------------------------------------------------------------------------
    # create treeSummarizedExperiment
    rowData(se) <- rowData(se)[, colnames(rowData(se)) != "nodeLab",
                               drop = FALSE]
    tse <- new("treeSummarizedExperiment", se,
               linkData = linkA, treeData = tree)

    return(tse)
}

