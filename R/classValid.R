#' valid treeSummarizedExperiment class
#'
#' \code{validTSE} is to valid treeSummarizedExperiment object.
#'
#'
#' @param obj A treeSummarizedExperiment object
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData
#' @importFrom S4Vectors metadata
#' @return TRUE or a character string.
#' @keywords internal

validTSE <- function(obj){

    # -------------------------------------------------------------------------
    # it must have table in assays
    if (length(assays(obj)) < 0){
        return("\n there is nothing in assays. \n")
    }

    # -------------------------------------------------------------------------
    #  Tree should be a phylo object
    if(!inherits(metadata(obj)$tree, "phylo")) {
        return("\n tree is not a phylo object")
    }

    # -------------------------------------------------------------------------
    # all leaf labels should be unique
    leafLab <- metadata(obj)$tree$tip.label
    isDup <- any(duplicated(leafLab))
    if(isDup){
        return("\n The leaf label of tree is not unique")
    }

    # -------------------------------------------------------------------------
    # if rownames exist, they should match with the node labels of the tree
    rName <- rownames(obj)
    labTree <- c(metadata(obj)$tree$tip.label, metadata(obj)$tree$node.label)
    if (!is.null(rName)) {
        isSame <- all(rName %in% labTree)
        if (!isSame) {
            difRN <- rName[!rName %in% labTree]
            stop("\n mismatch between assays data and tree: ",
                 head(difRN), " could not be found from tree. \n")
        }
    }

    # -------------------------------------------------------------------------
    # if rownames doesn't exist, two columns, nodeLab & rowID are required in
    # rowData.
    # if rowNum doesn't exist, the value in column nodeLab should match with
    # the node labels of tree.
    if (is.null(rName)) {
        # existance
        AA <- rowData(obj)$nodeLab
        BB <- rowData(obj)$rowID
        CC <- rowData(obj)$nodeNum
        isAA <- is.null(AA)
        isBB <- is.null(BB)
        if (any(isAA, isBB)) {
            stop("\n Two columns, nodeLab & rowID must be included in rowData
                 when rownames isn't provided.")
        } else {
            # match
            isCC <- is.null(CC)
            isSame <- all(AA %in% labTree)
            if (isCC & !isCC) {
                difAA <- setdiff(AA,labTree)
                stop("\n mismatch between nodeLab column (in rowData)
                     and tree: ", head(difAA), " could not be found from
                     tree. \n")
            }
            }
            }

    # -------------------------------------------------------------------------
    # Note : duplicated value in nodeLab column is allowed because, we might
    # have, for DS data, multiple rows correspondind to a same node.
    return(TRUE)
        }


