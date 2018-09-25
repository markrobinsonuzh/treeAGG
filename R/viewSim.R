#' visualize simulated scenario
#'
#' \code{viewSim} is to visualize the output from the function \code{simData}.
#'
#' @param obj The output from \code{simData}
#' @param layout The tree layout. Select one from 'rectangular', 'slanted',
#'   'fan', 'circular', 'radial', 'equal_angle' and 'daylight'. The default is
#'   "rectangular".
#' @param zoomScale A positive numeric value. If it is above one, branches with
#'   fold change equal to one (non-signal branch) will be zoomed in; If below
#'   one, they will be shrinked. Default is 0.05
#' @param legend.theme A list of arguments used for the theme in ggplot2
#' package (see \code{\link[ggplot2]{theme}} ) and starting with "legend."
#' @param tip.label TRUE or FALSE. Default is FALSE. If TRUE, the leaves with
#' fold change above or below 1 will be labelled.
#' @param legend.title The title of the legend. The default is "Abundance"
#'
#' @importFrom S4Vectors metadata
#' @importFrom ggtree %<+% geom_tiplab geom_point2
#' @export
#'
#' @return a figure
#'
viewSim <- function(obj, layout = "rectrangular", zoomScale = 1/20,
                    legend.theme = list(
                        legend.position = c(0.15, 0.6)),
                    tip.label = FALSE, legend.title = "Abundance"){

    md <- metadata(obj)
    # tree
    tree <- treeData(obj)

    # branch
    branch <- c(md$branch$A, md$branch$B)
    col.branch <- c( "go up" = "orange",
                     "go down" = "blue")

    # legend
    fc <- c(md$branch$"A_prop(%)", md$branch$"B_prop(%)")
    vv <- sum(diff(fc) > 0) + 1 # TRUE take 2; FALSE take 1
    nam <- switch(vv,
                  c("go down", "go up"),
                  c("go up", "go down"))
    legend.label <- list(col.branch = nam)
    legend.title <- c(branch = legend.title)

    # branch color
    colL <- col.branch[nam]

    # scenario
    sc <- md$scenario

    # zoomNode
    sNode <- findParallel(tree = tree, input = branch, label = FALSE)
    isT <- isLeaf(tree = tree, input = sNode)
    fNode <- sNode[!isT]

    # fold change
    d <- data.frame(node = transNode(tree = tree,
                                     input = names(md$FC)),
                    fc = md$FC)

    # figure
    if(sc == "S1") {
        fig <- treePlot(tree = tree, branch = branch,
                        col.branch = colL, zoomNode = fNode,
                        zoomScale = zoomScale, layout = layout,
                        legend = TRUE, legend.label = legend.label,
                        legend.theme = legend.theme,
                        legend.title = legend.title)

        if(tip.label) {
            fig <- fig %<+% d +
                geom_tiplab(aes(subset = (fc !=1)))
        }
    }
     if(sc == "S2") {
         fig <- treePlot(tree = tree, branch = branch,
                         col.branch = colL, col.point = "NA",
                         zoomNode = fNode,
                         zoomScale = zoomScale, layout = layout,
                         legend = TRUE,
                         legend.label = legend.label,
                         legend.theme = legend.theme,
                         legend.title = legend.title)

         fig <- fig %<+% d +
             geom_point2(aes(subset = (fc != 1),
                             size = fc))

         if(tip.label) {
             fig <- fig +
                 geom_tiplab(aes(subset = (fc !=1)))
         }
     }
    #
    if(sc == "S3") {
        difN <- names(md$FC)[md$FC != 1]
        difN <- transNode(tree = tree, input = difN)
        difS <- signalNode(node = difN, tree = tree)
        difL <- transNode(tree = tree, input = difS)
        colL2 <- ifelse(difL == md$branch$B,
                     colL[2], colL[1])
        ind <- match(colL2, col.branch)
        legend.label <- list(col.branch = names(col.branch)[ind],
                             col.other = "same")
        fig <- treePlot(tree = tree, branch = difS,
                        col.branch = colL2, zoomNode = fNode,
                        zoomScale = zoomScale, layout = layout,
                        legend = TRUE,
                        legend.label = legend.label,
                        legend.theme = legend.theme,
                        legend.title = legend.title)
        if(tip.label) {
            fig <- fig %<+% d +
                geom_tiplab(aes(subset = (fc !=1)))
        }
    }



    fig
}

#' find parallel nodes
#'
#' \code{findParallel} is to find target nodes that could build up a tree with
#' some specified nodes. The returned target nodes is the combination that has
#' the minimum number of nodes. In other words, a tree is cut in a way that
#' the number of branches is the minimum in all possibilities and the specified
#' branches are obtained. A branch is represented by its branch node. A leaf
#' node represents the edge connecting the leaf and its parent.
#'
#' @param tree A phylo object.
#' @param input A numeric or character vector. Node labels or node numbers.
#' @param label A logical value, TRUE or FALSE. If TRUE, node labels are
#' returned.
#'
#' @keywords internal
#' @return A vector of node labels or node numbers
#'
findParallel <- function(tree, input, label = FALSE){
    # find descendant leaves of input
    inT <- lapply(input,
                  FUN = function(x){
                      findOS(tree = tree,
                             ancestor = x,
                             only.Tip = TRUE)})
    inT <- unlist(inT)
    # find all leaves of tree
    allT <- setdiff(tree$edge[,2], tree$edge[, 1])

    # Leaves not included in input
    exT <- setdiff(allT, inT)

    # replace leaves with their ancestor branch node
    fT <- signalNode(tree = tree, node = exT, label = FALSE)

    # isT <- isLeaf(tree = tree, input = fT)
    # fn <- fT[!isT]
    return(fT)
}

#' To test whether specified nodes are leaf nodes
#'
#' \code{isLeaf} is to test wheter some specified nodes are leaf nodes of a tree.
#'
#' @param tree A phylo object.
#' @param input A numeric or character vector. Node labels or node numbers.
#'
#' @return a logical vector with the same length as input
#' @keywords internal
isLeaf <- function(tree, input){
    # leaves
    tip <- setdiff(tree$edge[,2], tree$edge[, 1])
    # input
    if (inherits(input, "character")) {
        input <- transNode(tree, input = input)
    } else {
        input <- input
    }
    # is it a tip
    isLeaf <- input %in% tip
    return(isLeaf)
}
