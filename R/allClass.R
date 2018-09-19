
# ======================
### class
# ======================

#  -----------------------------------------------------------------------------
#' int_rowData: A virtual class
#' -----------------------------------------------------------------------------
#' A virtual class to be assigned to a column of \code{rowData} in a
#' \code{treeSummarizedExperiment} object. The columns with this class could be
#' optionally exported or hidden when using \code{rowData}.
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
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} class. It is
#' designed to store rectangular data like a \code{matrix} for entities (e.g.,
#' microbes or cell types), and the hiearchical structure information of
#' entities. The annotations on the rows and columns of the rectangular data are
#' strored in \code{rowData} and \code{colData}, respectively. Each row of the
#' rectangular data could be mapped to a leaf node of the tree via the rownames
#' of \strong{rowData} or a column named as \strong{nodeLab} in
#' \strong{rowData}. For example, we could include abundance count of microbes
#' collected from different samples in the \strong{assays} and the phylogenetic
#' tree of microbes in the \strong{metadata}. Each row of matrix-like data in
#' \strong{assays} corresponds to one microbial specie, whose information is
#' given in \strong{rowData}. We could find the label of its corresponding leaf
#' node in column \strong{nodeLab}, or in the rownames of \strong{rowData}. Each
#' column of matrix-like data is a sample, whose information is stored in
#' \strong{colData}. The class \strong{leafSummarizedExperiment} has four slots
#' \strong{assays}, \strong{rowData} \strong{colData} and \strong{metadata} as
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} class. More details
#' about these four slots could be found in
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#'
#'
#' \strong{leafSummarizedExperiment} has more restrictions in data structure
#' than \strong{SummarizedExperiment}.
#' In \strong{leafSummarizedExperiment} class, we have restrictions as below.
#' \itemize{
#' \item It's not allowed to have empty \strong{assays}.
#'
#' \item A phylo object named as \strong{tree} should be contained in
#' \strong{metadata}. The labels of tree leaf nodes are unique.
#'
#' \item Each rows of matrix-like data in \strong{assays} could be mapped to a
#' leaf node of the tree. The label of leaf node is provided in a column
#' \strong{nodeLab} in \strong{rowData}, or provides as the rownames of
#' \strong{rowData}. If both provided, the information in column \strong{nodeLab} is
#' used. }
#'
#' @section Constructor:
#' See \code{\link{leafSummarizedExperiment-constructor}} for constructor
#' functions.
#'
#' @section Accessor:
#' See \code{\link{leafSummarizedExperiment-accessor}} for accessor functions.
#'
#' @importFrom methods setClass
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @name leafSummarizedExperiment-class
#' @exportClass leafSummarizedExperiment
#' @include classValid.R
.lse <- setClass("leafSummarizedExperiment",
                 contains = "SummarizedExperiment",
                 validity = checkLSE)

#-------------------------------------------------------------------------------
# treeSummarizedExperiment
#-------------------------------------------------------------------------------
#' An S4 class treeSummarizedExperiment
#'
#' The class \strong{treeSummarizedExperiment} is an extension class of standard
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} class. It has six
#' slots. Four of them are traditional slots from
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} class:
#' \strong{assays}, \strong{rowData} \strong{colData} and \strong{metadata}. The
#' other two slots are \strong{linkData} and \strong{tree}. The class
#' \strong{treeSummarizedExperiment} is designed to store rectangular data for
#' entities (e.g., microbes or cell types) (\strong{assays}), information about
#' the hiearchical structure of entities (\strong{tree}), and information
#' about the mapping between the rectangular data and the tree
#' (\strong{linkData}).
#'
#' @slot linkData A \code{\link[S4Vectors]{DataFrame}} object. It gives map
#'   information between the rows of rectangular data and the nodes of tree.
#'   The row dimension is the same to that of \code{assays} data
#'   \itemize{
#'   \item \strong{nodeLab} the node labels on the tree.
#'   \item \strong{nodeLab_alias} a alias of column \code{nodeLab}. It is
#'   created only when there are missing value or duplicated value in column
#'   \code{nodeLab}. Commonly, it is a prefix "Node_" follow by the value in
#'   column \code{nodeNum}.
#'   \item \strong{nodeNum} the node numbers on the tree.
#'   \item \strong{isLeaf} is it a leaf node?
#'   \item \strong{rowID} the row number in \code{assays}.
#'   }
#' @slot treeData A phylo object. It gives information about the hiearchical
#'   structure of entities.
#' @slot ... See \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'   for more details about slots inherited from \code{SummarizedExperiment}
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




# ======================
### builder
# ======================
#' constructor for leafSummarizedExperiment object
#'
#' \code{leafSummarizedExperiment} is to create a
#' \strong{leafSummarizedExperiment} object.
#'
#' @param tree A phylo object.
#' @param ... see arguments in
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' \itemize{
#' \item{\strong{assays}:}{ a list of matrix-like elements, or a matrix-like
#' object. All elements of the list must have the same dimensions, and
#' dimension names (if present) must be consistent across elements and with
#' the row names of \strong{rowData} and \strong{colData}}.
#'
#' \item{\strong{rowData}:}{ A \code{\link[S4Vectors]{DataFrame-class}} object
#' describing the rows. Row names, if present, become the row names of the
#' \strong{leafSummarizedExperiment}. The number of rows of the
#' \code{\link[S4Vectors]{DataFrame-class}} must equal the number of rows of the
#' matrices in \strong{assays}. Note: a column named \strong{nodeLab} or
#' rownames, which provides the labels of corresponding nodes on the tree for
#' each row of the matrices in \strong{assays}, is required.}
#'
#' \item{\strong{colData}:}{ An optional
#' \code{\link[S4Vectors]{DataFrame-class}} describing the samples. Row names,
#' if present, become the column names of the
#' \strong{leafSummarizedExperiment}.}
#'
#' \item{\strong{metadata}:}{ An optional list of arbitrary content describing
#' the overall experiment.}
#' }
#'
#' @importFrom utils head
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData
#'   "rowData<-" "metadata<-"
#' @importFrom S4Vectors metadata
#' @importFrom methods as
#' @export
#' @return A leafSummarizedExperiment object
#' @seealso \code{\link{leafSummarizedExperiment-class}}
#'   \code{\link{SummarizedExperiment-class}}
#' @name leafSummarizedExperiment-constructor
#' @author Ruizhu HUANG
#' @examples
#' # tree
#' data("tinyTree")
#'
#' # assays
#' count <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count) <- c(tinyTree$tip.label)
#' colnames(count) <- paste("C_", 1:10, sep = "_")
#'
#' # colData
#' sampC <- data.frame(condition = rep(c("control", "trt"), each = 5),
#' gender = sample(x = 1:2, size = 10, replace = TRUE))
#' rownames(sampC) <- colnames(count)
#'
#' lse <- leafSummarizedExperiment(tree = tinyTree, assays = list(count),
#' colData = sampC)
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
            cat(sum(isOut), "rows are removed from tables of *assays*. \n",
                nodeLab[isOut]," cannot match to any label of tree leaf node. \n")}
        se <- se[isIn, ] }
    # use rownames
    if (is.null(nodeLab)) {
        rowNam <- rownames(se)
        if (is.null(rowNam)) {
            stop("Either rownames or a nodeLab column in rowData should be
                 provided \n.")
        }
        isIn <- rowNam %in% tipLab
        isOut <- !isIn
        if (sum(isOut) > 0) {
            cat(sum(isOut), "rows are removed from tables of *assays*. \n",
                rowNam[isOut], " cannot match with any label of tree leaf node. \n")}
        se <- se[isIn, ]
        }

    # output as leafSummarizedExperiment
    #lse <- as(se, "leafSummarizedExperiment")
    lse <- .lse(se)
    return(lse)
    }

#' A Constructor for treeSummarizedExperiment
#'
#' \code{treeSummarizedExperiment} constructs a treeSummarizedExperiment.
#'
#' @param tree A phylo object
#' @param linkData A data frame or \code{\link[S4Vectors]{DataFrame}}. It has
#'   the same number of rows as the matrix-like element in \code{assays.}. The
#'   row order between \code{linkData} and the matrix-like element in
#'   \code{assays.} is the same. It includes columns as below.
#'   \itemize{
#'   \item \strong{nodeLab} the labels of nodes
#'   \item \strong{nodeLab_alias} a alias of column \code{nodeLab}. It is
#'   created only when there are missing value or duplicated value in column
#'   \code{nodeLab}. Commonly, it is a prefix "Node_" follow by the value in
#'   column \code{nodeNum}.
#'   \item \strong{nodeNum} the numbers of nodes
#'   \item \strong{isLeaf} is it a leaf node
#'   \item \strong{rowID} the row number in the matrix-like element in
#'   \code{assays} }
#' @param ... see arguments in
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' \itemize{
#' \item{\strong{assays}:}{ a list of matrix-like elements, or a matrix-like
#' object. All elements of the list must have the same dimensions, and
#' dimension names (if present) must be consistent across elements and with
#' the row names of \strong{rowData} and \strong{colData}}.
#'
#' \item{\strong{rowData}:}{ A \code{\link[S4Vectors]{DataFrame-class}} object
#' describing the rows. Row names, if present, become the row names of the
#' \strong{leafSummarizedExperiment}. The number of rows of the
#' \code{\link[S4Vectors]{DataFrame-class}} must equal the number of rows of the
#' matrices in \strong{assays}. }
#'
#' \item{\strong{colData}:}{ An optional
#' \code{\link[S4Vectors]{DataFrame-class}} describing the samples. Row names,
#' if present, become the column names of the
#' \strong{leafSummarizedExperiment}.}
#'
#' \item{\strong{metadata}:}{ An optional list of arbitrary content describing
#' the overall experiment.}
#' }
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
#' @export
#' @return a treeSummarizedExperiment object
#' @name treeSummarizedExperiment-constructor
#' @author Ruizhu HUANG
#' @seealso \code{\link{treeSummarizedExperiment}}
#'   \code{\link{treeSummarizedExperiment-accessor}}
#'   \code{\link{leafSummarizedExperiment}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
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
                cat(sum(isOut), "rows are removed. They cannot be assigned to
                    any tree nodes. \n")}
            se <- se[isIn, ]
            newLab <- rowData(se)$nodeLab
            }

        # (2) if the nodeLab column doesn't exist, rownames should match with the
        # labels of tree leaves.
        if (is.null(nodeLab)) {
            rowNam <- rownames(se)

            # if neither nodeLab nor rownames are provided, return error.
            if (is.null(rowNam)) {
                stop("Either rownames or a nodeLab column in rowData should be
                     provided \n.")}
            # keep only rows that could be assigned to the nodes of tree
            isIn <- rowNam %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                cat(sum(isOut), "rows are removed. They cannot be assigned to
                    any tree nodes. \n")}
            se <- se[isIn, ]
            newLab <- rownames(se)
            }

        # create linkD
        linkD <- DataFrame(nodeLab = newLab,
                           nodeNum = transNode(tree = tree,
                                               input = newLab),
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
    rowData(se) <- rowData(se)[, colnames(rowData(se)) != "nodeLab"]
    tse <- new("treeSummarizedExperiment", se,
               linkData = linkD, treeData = tree)

    return(tse)



}
