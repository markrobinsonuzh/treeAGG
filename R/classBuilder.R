

#' create a treeSummarizedExperiment object
#'
#' \code{treeSummarizedExperiment} is to create a object with
#' \strong{treeSummarizedExperiment} class. The class
#' \strong{treeSummarizedExperiment} is an extension to
#' \strong{SummarizedExperiment}.
#'
#' @param tree A phylo object.
#' @param ... see arguments for
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' \itemize{
#' \item{\strong{assays}:}{ a list of matrix-like elements, or a matrix-like
#' object. All elements of the list must have the same dimensions, and
#' dimension names (if present) must be consistent across elements and with
#' the row names of \strong{rowData} and \strong{colData}}.
#' \item{\strong{rowData}:}{ A \code{\link[S4Vectors]{DataFrame-class}}
#' describing the rows. Row names, if present, become the row names of the
#' \strong{treeSummarizedExperiment}. The number of rows of the DataFrame must
#' equal the number of rows of the matrices in \strong{assays}. Note: a column
#' named \strong{nodeLab} is required to provide the corresponding nodes on the
#' tree for each row of the matrices in \strong{assays}.}
#' \item{\strong{colData}:}{ An optional DataFrame describing the samples.
#' Row names, if present, become the column names of the
#' \strong{treeSummarizedExperiment}.}
#' \item{\strong{metadata}:}{ An optional list of arbitrary content describing the overall
#'  experiment.}}
#' @importFrom utils head
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData
#'   "rowData<-" "metadata<-"
#' @importFrom methods as
#' @export
#' @return A treeSummarizedExperiment object
#' @seealso \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'
#' @examples
#' # tree
#' data("tinyTree")
#'
#' # assays
#' count <- matrix(rpois(190, 50), nrow = 19)
#' rownames(count) <- c(tinyTree$tip.label, tinyTree$node.label)
#' colnames(count) <- paste("C_", 1:10, sep = "_")
#'
#' # colData
#' sampC <- data.frame(condition = rep(c("control", "trt"), each = 5),
#' gender = sample(x = 1:2, size = 10, replace = TRUE))
#' rownames(sampC) <- colnames(count)
#'
#' tse <- treeSummarizedExperiment(tree = tinyTree, assays = list(count),
#' colData = sampC)
#'
treeSummarizedExperiment <- function(tree, ...){

    # -------------------------------------------------------------------------
    # tree is a phylo object
    if(!inherits(tree, "phylo")) {
        stop(tree, ": A phylo object is required", "\n")
    }
    treeLab <- c(tree$tip.label, tree$node.label)
    # -------------------------------------------------------------------------
    # create a SummarizedExperiment object
    obj <- SummarizedExperiment(...)

    # -------------------------------------------------------------------------
    # generate or update rowData
    asy <- assays(obj)
    rd <- rowData(obj)
    nrd <- ncol(rd)
    rn <- rownames(obj)
    isLeaf <- rn %in% tree$tip.label

    if(length(asy) > 0){
        isRN <- is.null(rn)
        isNL <- is.null(rd$nodeLab)

        # match table with tree:
        # 1. if only rownames exists, use rownames
        if((!isRN) & isNL){
            isSame <- all(rn %in% treeLab)
            if(!isSame){
                diffLab <- rn[! rn %in% treeLab]
                stop("\n Mismatch between table and tree, rownames: ",
                     head(diffLab), " can't be found on tree. \n")
            }
            # add columns into rowData: nodeLab, isLeaf, nodeNum
            rowData(obj)$nodeLab <- rn
            rowData(obj)$nodeNum <- transNode(tree = tree,
                                              input = rn)
            rowData(obj)$isLeaf <- isLeaf
            rowData(obj)$rowID <- seq_along(rn)
        }

        # 2. if nodeLab column (in rowData) exists, use nodeLab
           # (1) check values in nodeLab
           #     - If there are NA values, check the existance of column
           #       nodeNum. Return error if nodeNum isn't provided.
           #     - If there are not NA values, check whether all values in
           #       nodeLab exist in the tree. Return error if any doesn't exist.
        if(!isNL){
            nodeLab <- rd$nodeLab
            if(!inherits(nodeLab, "character")){
                stop("\n The nodeLab column in rowData should be character. \n")
            }

            # if there are NA values in nodeLab, check the existance of nodeNum
            isNA <- any(is.na(nodeLab))
            if (isNA) {
                isNN <- is.null(rd$nodeNum)
                if (isNN) {
                    stop("\n NA value in nodeLab; Column nodeNum is required to
                         to match table rows and tree nodes. \n")
                }
            } else {
            # if not NA values, check whether all values exist in tree
            isSame <- all(nodeLab %in% treeLab)
            if(!isSame){
                diffLab <- nodeLab[! nodeLab %in% treeLab]
                stop("\n Mismatch between nodeLab column and tree node
                     labels: ", head(diffLab), " can't be found on tree. \n")
            }

            # add columns into rowData: isLeaf, nodeNum
            if(is.null(rd$nodeNum)){
                rowData(obj)$nodeNum <- transNode(tree = tree,
                                                  input = rd$nodeLab)
            }
            if(is.null(rd$isLeaf)){
                rowData(obj)$isLeaf <- isLeaf
            }
            if(is.null(rd$rowNum)){
                rowData(obj)$rowID <- seq_along(nodeLab)
            }

        }
        }
        # 3. if neither exist, return error.
        if(isRN & isNL){stop("\n The information to match table with tree is
                             not given \n")}

        # add tree in metadata
        metadata(obj)$tree <- tree
        new_obj <- as(obj, "treeSummarizedExperiment")

        }else{
            stop("\n assays data is missing... \n")
        }
    #validObject(new_obj)
    return( new_obj )

}


