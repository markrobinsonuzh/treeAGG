#' visualize the phylogenetic tree
#'
#' \code{treePlot} is to visualize a phylogenetic tree.
#'
#' @param tree a phylo object
#' @param branch a vector of node numbers or node labels to specify the branches to be colored. Each branch is represented by its root node(node numbers or labels). A leaf node reprents the edge connecting the leaf and its parent.
#' @param col.branch a vector of colors. Its length should be one or equals to the length of \strong{branch}. If the lengths of col.branch and branch are the same, their corresponding orders are matched. The default is blue.
#' @param col.other a color for the branches other than those specified in branch
#' @param point a vector of node numbers or node labels to specify the point locations in the tree
#' @param col.point the color for the \strong{point}.
#' @param size.point the size for the \strong{point}.
#' @param zoomNode a vector of nodes to be zoomed in. If default (null), the tree is not zoomed in.
#' @param zoomLevel a numeric vector. Its length is equal to 1 or equal to the length of zoomNode. If default (null), a leaf is zoomed in its direct parent level and a internal node is zoomed in its own level.
#'
#' @param zoomScale a scale to zoom in. If default (null), tree is not zoomed in.
#' @param legend.title a vector to specify the title of the legend. It must be named with "branch" and "point" to match with the argument \strong{branch} and \strong{point}.
#' @param size.line.legend the line size shown in the legend for \strong{branch}
#' @param size.point.legend the point size shown in the legend for \strong{point}.
#' @param ... see also \code{\link[ggtree]{ggtree}}
#'
#' @details treePlot is created based on the \pkg{ggtree} and \pkg{ggplot2}. So other geoms from these two packages could be combined with \code{treePlot} to add geoms in the figure created by \code{treePlot}.
#'
#' @export
#' @import ggplot2 ggtree
#' @return a tree plot
#'
#' @examples
#'
#' data(exTree)
#' ggtree(bigTree, layout = "circular", branch.length = "none")
#'
#' # If we want to color two branches with root node 1000 and 1400
#' p <- treePlot(tree = bigTree, branch = c(1000, 1400))
#' p
#'
#' # To use different color
#' p <- treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("red", "blue"), col.other = c("other" = "grey20"))
#' p
#'
#' # To add legend, we could do it in this way
#' # the name of col.other would be used as the label in the legend
#' p + ggplot2::theme(legend.position = "right")
#'
#' # To change the legend title
#' p <- treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("red", "blue"), col.other = c("other" = "grey20"),
#' legend.title = c("branch" = "Truth"))
#'
#' p + ggplot2::theme(legend.position = "right")
#'
#' # To zoom in the branch with root node 1400 (by scale 4)
#' p <- treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("red", "blue"), col.other = c("other" = "grey20"),
#' legend.title = c("branch" = "Truth"), zoomNode = 1400,
#' zoomScale = 8)
#'
#' p + ggplot2::theme(legend.position = "right")
#'
#' # If its sister branch is also interested, we could zoom in their parent level by zoomLevel = 1
#'
#' p <- treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("red", "blue"), col.other = c("other" = "grey20"),
#' legend.title = c("branch" = "Truth"), zoomNode = 1400,
#' zoomScale = 8, zoomLevel = 1)
#'
#' p + ggplot2::theme(legend.position = "right")
#'
#' # If we want to add some points at some nodes of the tree
#' p <- treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("red", "blue"), col.other = c("other" = "grey20"),
#' legend.title = c("branch" = "Truth"), zoomNode = 1400,
#' zoomScale = 8, zoomLevel = 1, point = c(50, 100, 500, 1200, 1500))
#'
#' # To change the legend for the point
#' p <- treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("coral", "blue"), col.other = c("other" = "darkgrey"),
#' legend.title = c("branch" = "Truth", "point" = "Estimate"),
#'  zoomNode = 1400,
#' zoomScale = 8, zoomLevel = 1, point = c(50, 100, 500, 1200, 1500),
#' col.point = c("found" = "violet"))
#'
#' p + ggplot2::theme(legend.position = "bottom",
#' legend.text = element_text(size= 10),
#' legend.key.size = unit(4,"cm"),
#' legend.key.height = unit(0.4,"cm"),
#' legend.key.width = unit(0.5, "cm"),
#' legend.title = element_text(size = 15))

treePlot <- function(tree, branch = NULL,
                     col.branch = "blue",
                     col.other = "grey",
                     point = NULL,
                     col.point = "orange",
                     size.point = 2,
                     zoomNode = NULL,
                     zoomLevel = NULL,
                     zoomScale = NULL,
                     legend.title = c("branch"= "Title_branch",
                                      "point" = "Title_point"),
                     size.line.legend = 2,
                     size.point.legend = 2,...){
  # check tree
  if (!inherits(tree, "phylo")) {
    stop("tree: should be a phylo object")
  }

  p <- ggtree(tree)
  d <- p$data[, "node", drop = FALSE]

  # ================== points ===============================
  # add points on the tree
  if(!is.null(point)){

    # transformation from node label to node number
    if (inherits(point, "character")) {
      point <- tx_node(tree = tree, input = point)
    } else {
      point <- point
    }

    if (is.null(names(col.point))) {
      names(col.point) <- "col.point"
    } else {
      col.point <- col.point
    }

    # color
    colG <- col.point

    # create a data frame to store the information for points
    d$Estimate <- ifelse(d$node %in% point, names(col.point),
                         "NO_Found")
    d$show <- ifelse(d$node %in% point, TRUE, FALSE)

    # legend title
    if(is.na(legend.title["point"])){
      legend.title["point"] <- ""
    }

    p1 <- ggtree(tree, ...) %<+% d +
      geom_point2(aes(subset = show, color = Estimate,
                      size = Estimate)) +
      scale_size_manual(values = size.point,
                        guide = guide_legend(
                          title = legend.title["point"],
                          override.aes = list(
                            shape = 16, color = col.point,
                            size = size.point.legend)))

  }else{
    # color
    colG <- c()
    p1 <- ggtree(tree, ...)

  }

  # ================== edges ===============================
  if(!is.null(branch)){

    # transformation from node label to node number
    if (inherits(branch, "character")) {
      branch <- tx_node(tree = tree, input = branch)
    } else {
      branch <- branch
    }

    # check to see if col.branch has propriate length
    if(length(col.branch) != length(branch) &
       length(col.branch) != 1) {
      stop("The length of col.branch should be one
           or equal to the length of branch")
    }

    # organize the color used in the figure
    if(is.null(names(col.branch))) {
      if (length(col.branch) != length(branch)) {
        names(col.branch) <- "col.branch"
      } else { names(col.branch) <- branch}
    } else {
      col.branch <- col.branch }

    if(is.null(names(col.other))) {
      names(col.other) <- "col.other"
    } else {
      col.other <- col.other}

    colG <- c(colG, col.branch, col.other)


    # The edges selected to be colored
    desList <- lapply(branch,
                      FUN = function(x) {
                        findOS(ancestor = x, tree = tree,
                               only.Tip = FALSE, self.include = TRUE)})
    names(desList) <- branch
    des_len <- unlist(lapply(desList, length))
    des.1 <- desList[order(des_len, decreasing = TRUE)]
    group_var <- lapply(seq_along(des.1), FUN = function(x) {
      rep(names(des.1)[x], length(des.1[[x]]))})
    df.1 <- data.frame(node = unlist(des.1, use.names = FALSE),
                       group = unlist(group_var),
                       stringsAsFactors = FALSE)

    if (length(col.branch) == 1) {
      df.1$group <- names(col.branch)
    }

    # create a data frame to include the selection information for edges
    d$Truth <- rep(names(col.other), nrow(d))
    d$Truth[match(df.1$node, d$node)] <- df.1$group

    # integrate with tree
    if(is.na(legend.title["branch"])){
      legend.title["branch"] <- ""
    }
    p1 <- p1 %<+% d + aes(colour = Truth) +
      scale_color_manual(values = colG,
                         labels = ifelse(! names(colG) %in% names(col.point),
                                         names(colG), ""),
                         guide = guide_legend(
                           title = legend.title["branch"],
                           override.aes =list(
                             color = colG, linetype = ifelse(
                               names(colG) %in% names(col.point),
                               "blank", "solid"),
                             shape = rep(NA, length(colG)),
                             size = size.line.legend)))


    }else{
      p1 <- p1 +
        scale_color_manual(values = colG,
                           labels = ifelse(!names(colG) %in% names(col.point),
                                           names(colG), ""),
                           guide = FALSE)
    }

  # ================ ZOOM IN ====================
  #
  if(!is.null(zoomNode)){

    if (inherits(point, "character")) {
      zoomNode <- tx_node(tree = tree, input = zoomNode)
    } else {
      zoomNode <- zoomNode
    }

    zList <- lapply(zoomNode,
                    FUN = function(x) {
                      findOS(ancestor = x, tree = tree,
                             only.Tip = FALSE, self.include = TRUE)})
    names(zList) <- zoomNode
    z_len <- unlist(lapply(zList, length))

    # define zoomLevel
    if (is.null(zoomLevel)) {
      zoomLevel <- ifelse(z_len > 1, 0, 1)
    } else {
      if (length(zoomLevel) == 1) {
        zoomLevel <- rep(zoomLevel, length(zoomNode))
      } else {
        zoomLevel <- zoomLevel
      }
    }
    names(zoomLevel) <- zoomNode

    # define zoomScale
    if (is.null(zoomScale)) {
      zoomScale <- rep(1, length(zoomNode))
    } else {
      zoomScale <- rep(zoomScale, length(zoomNode))
    }

    # the nodes to be zoomed in
    nodZ <- findAncestor(tree = tree, node = zoomNode, level = zoomLevel)
    names(nodZ) <- names(zoomScale) <- zoomNode
    # remove nodes which are the descendants of the others
    nodZW <- rmMember(node = nodZ, tree = tree)
    zoomScale[!names(zoomScale) %in% names(nodZW)] <- 1

    # zoom the selected nodes
    i <- 1
    repeat {
      p1 <- p1 %>% scaleClade(nodZ[i], scale = zoomScale[i])
      i <- i + 1
      if (i > length(nodZ)) {
        break
      }
    }


    lim <- c(min(p1$data$y), max(p1$data$y))

    ## ggtree function set ylim when layout is circular or radical this would lead to
    ## issue, like points not displayed when zoom in some branches
    suppressMessages(p1 <- p1 + scale_y_continuous(limits = lim))

  }

  return(p1)
}


