
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
treeSummarizedExperiment <- function(tree = NULL, linkData = NULL,
                                     ...){

    # -------------------------------------------------------------------------
    # create a SummarizedExperiment object
    se <- SummarizedExperiment(...)

    # -------------------------------------------------------------------------
    # tree is a phylo object
    if(!inherits(tree, "phylo")) {
        stop(tree, ": A phylo object is required", "\n")
    }

    # -------------------------------------------------------------------------
    # The labels of tree nodes should be unique
    treeLab <- c(tree$tip.label, tree$node.label)
    tipLab <- tree$tip.label
    isDp <- duplicated(tipLab)
    anyDp <- any(isDp)
    if (anyDp) {
        stop("\n Can not distinguish nodes with the same label: ",
             head(tipLab[isDp])," \n")
    }

    # -------------------------------------------------------------------------
    ### create link data
    if (is.null(linkData)) {

        # (1) if the nodeLab column exist, they should match with the labels of
        # tree nodes.
        nodeLab <- rowData(se)$nodeLab
        if (!is.null(nodeLab)) {
            # keep only rows that could be assigned to the nodes of tree
            isIn <- nodeLab %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                message(sum(isOut), "rows are removed. They can't be matched to
                    any node of the tree. \n")}
            se <- se[isIn, ]
            newLab <- rowData(se)$nodeLab
            }

        # (2) if the nodeLab column doesn't exist, rownames should match with
        # the labels of tree leaves.
        if (is.null(nodeLab)) {
            rowNam <- rownames(se)

            # if neither nodeLab nor rownames are provided, return error.
            if (is.null(rowNam)) {
                stop("Either a nodeLab column or row names should be
                     provided for row data \n.")}
            # keep only rows that could be assigned to the nodes of tree
            isIn <- rowNam %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                message(sum(isOut), " rows are removed. They cannot be maped to
                    any node of the tree. \n")}
            se <- se[isIn, ]
            newLab <- rownames(se)
            }

        # create linkD
        linkD <- DataFrame(nodeLab = newLab,
                           nodeNum = transNode(tree = tree,
                                               input = newLab,
                                               message = FALSE),
                           isLeaf = newLab %in% tree$tip.label,
                           rowID = seq_len(length(newLab)))

        # create column nodeLab_alias, if there are duplicated value in the
        # nodeLab column
        if (any(duplicated(linkD$nodeLab))) {
            linkD$nodeLab_alias <- paste("Node_", linkD$nodeNum, sep = "")
        }

        } else {
            # if linkData is provided, then use it as linkData
            linkD <- linkData

            # create column nodeLab_alias, if there are duplicated value in the
            # nodeLab column
            if (any(duplicated(linkD$nodeLab))) {
                linkD$nodeLab_alias <- paste("Node_", linkD$nodeNum, sep = "")
            }
        }

    # -------------------------------------------------------------------------
    # create treeSummarizedExperiment
    rowData(se) <- rowData(se)[, colnames(rowData(se)) != "nodeLab",
                               drop = FALSE]
    tse <- new("treeSummarizedExperiment", se,
               linkData = linkD, treeData = tree)

    return(tse)



}
