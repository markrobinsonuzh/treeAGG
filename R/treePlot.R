#' Visualize the phylogenetic tree
#'
#' \code{treePlot} visualizes a phylogenetic tree.
#'
#' @param tree A phylo object
#' @param branch A vector of node numbers labels to specify the branches to be
#'   colored. Each branch is represented by its branch node. A leaf node
#'   reprents the edge connecting the leaf and its parent.
#' @param col.branch A vector of colors. Its length should be one or equals to
#'   the length of \strong{branch}. If \strong{col.branch} has the same length
#'   as \strong{branch}, the branches are colored correspondingly with the
#'   \strong{col.branch}. The default is blue.
#' @param col.other A color for the branches other than those specified in
#' \strong{branch}
#' @param point A vector of node numbers or node labels to specify the
#' locations to add points in the tree
#' @param col.point A color for the \strong{point}. It has length equal to one.
#' @param size.point The size for the \strong{point}. It has length equal to
#'   one.
#' @param zoomNode A vector of nodes to be zoomed in. If default (NULL), the
#'   tree is not zoomed in.
#' @param zoomLevel A numeric vector. Its length is equal to 1 or equal to the
#'   length of \strong{zoomNode}. If default (NULL), a leaf is zoomed in its
#'   direct parent level and an internal node is zoomed in its own level.
#' @param zoomScale A numeric vector. Its length is equal to one or equal to the
#'   length of \strong{zoomNode}. If \strong{zoomScale} has the same length as
#'   \strong{zoomNode}, the branches are zoomed in with different scales
#'   corresponding to the value of \strong{zoomScale}. If default (NULL), tree
#'   is not zoomed in.
#' @param legend TRUE or FALSE. Default is FALSE. If TRUE, the legend is
#'   created.
#' @param legend.theme A list of arguments used for the theme in ggplot2 package
#'   (see \code{\link[ggplot2]{theme}} ) and starting with "legend."
#' @param legend.title A vector to specify the title of the legend. It must be
#'   named with "branch" and "point" to match with the argument \strong{branch}
#'   and \strong{point}.
#' @param legend.label A list with three members: "col.branch", "col.other", and
#'   "col.point". The elements order in each member matches with the
#'   corresponding argument \strong{col.branch}, \strong{col.other} and
#'   \strong{col.point}, and will display in the legend.
#' @param size.line.legend The line size shown in the legend for \strong{branch}
#' @param size.point.legend The point size shown in the legend for
#' \strong{point}.
#' @param ... see also \code{\link[ggtree]{ggtree}}
#'
#' @details treePlot is created based on the \pkg{ggtree} and \pkg{ggplot2}. So
#'   other geoms from these two packages could be combined with \code{treePlot}
#'   to add geoms in the figure created by \code{treePlot}.
#'
#' @import ggplot2
#' @import ggtree
#' @export
#' @return A tree plot
#' @author Ruizhu Huang
#' @examples
#'
#' data(bigTree)
#'
#' # If we want to color two branches with branch node 1000 and 1400
#' treePlot(tree = bigTree, branch = c(1000, 1400))
#'
#'
#' # use col.branch and col.other to specify colors
#' treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("salmon", "blue"), col.other = "grey40")
#'
#' # add legend to the colored branches
#' treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("salmon", "blue"), col.other = "grey40",
#' legend = TRUE, legend.label = list(col.branch = c("up", "down")))
#'
#' # change legend title
#' p <- treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("salmon", "blue"), col.other = "grey40",
#' legend = TRUE,
#' legend.label = list(col.branch = c("Go up", "Go down")),
#' legend.title = c("branch" = "Abundance"))
#'
#' # change legend position (combine with ggplot2 package)
#' library(ggplot2)
#'  p + ggplot2::theme(legend.position = "bottom")
#'
#' # change legend position use legend.theme
#' treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("salmon", "blue"), col.other = "grey40",
#' legend = TRUE,
#' legend.label = list(col.branch = c("Go up", "Go down")),
#' legend.title = c("branch" = "Truth"),
#' legend.theme = list(legend.position = "bottom"))
#'
#'
#' # add points
#' treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("salmon", "blue"), col.other = "grey40",
#' legend = TRUE,
#' legend.label = list(col.branch = c("Go up", "Go down")),
#' legend.title = c("branch" = "Truth"),
#' legend.theme = list(legend.position = "bottom"),
#' point = c(500, 5, 10))
#'
#'
#'# add points label in legend
#' treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("salmon", "blue"), col.other = "grey40",
#' legend = TRUE,
#' legend.label = list(col.branch = c("Go up", "Go down"),
#' col.point = "Found"),
#' legend.title = c("branch" = "Truth", "point"= "Estimate"),
#' legend.theme = list(legend.position = "bottom"),
#' point = c(500, 5, 10))
#'
#'
#'# add points label in legend
#' treePlot(tree = bigTree, branch = c(1000, 1400),
#' col.branch = c("salmon", "blue"), col.other = "grey40",
#' legend = TRUE,
#' legend.label = list(col.branch = c("Go up", "Go down"),
#' col.point = "Found", col.other = "Same"),
#' legend.title = c("branch" = "Truth", "point"= "Estimate"),
#' legend.theme = list(legend.position = "bottom"),
#' point = c(500, 5, 10))
#'
treePlot <- function(tree,
                     branch = NULL,
                     col.branch = "blue",
                     col.other = "grey",
                     point = NULL,
                     col.point = "orange",
                     size.point = 2,
                     zoomNode = NULL,
                     zoomLevel = NULL,
                     zoomScale = 8,
                     legend = FALSE,
                     legend.theme = NULL,
                     legend.title = c(
                         "point" = "Title_point",
                         "branch" = "Title_branch"),
                     legend.label = NULL,
                     size.line.legend = 2,
                     size.point.legend = 3, ...) {

    # check tree
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }
    p <- ggtree(tree, ...)

    # color branch
    if (!is.null(branch)) {
        p <-  addBranch(tree = tree, branch = branch,
                        col.branch = col.branch,
                        col.other = col.other,
                        addTo = p)
    } else {
        p <- p
    }

    # add points
    if (!is.null(point)) {
        p <- addPoint(tree = tree, point = point,
                      col.point = col.point, addTo = p)
    } else {
        p <- p
    }

    # customize the size scale for added points
    if (!is.null(point)) {
        p <- p + sizeScale(size.point = size.point,
                           legend.label = legend.label,
                           legend.title = legend.title["point"],
                           col.point = col.point,
                           size.point.legend = size.point.legend,
                           legend = legend)
    } else {
        p <- p
    }

    # customize the color
    if (!is.null(branch)) {
        p <- p +
            colScale(branch = branch,
                     point = point,
                     col.branch = col.branch,
                     col.other = col.other,
                     col.point = col.point,
                     legend.label = legend.label,
                     legend.title = legend.title,
                     size.line.legend = size.line.legend,
                     legend = legend )
    } else {
        p <- p
    }

    # zoom in selected branches
    if (!is.null(zoomNode)) {
        p <- addZoom(tree = tree, zoomNode = zoomNode,
                     zoomLevel = zoomLevel, zoomScale = zoomScale,
                     addTo = p)
    } else {
        p <- p
    }

    # add legend
    if (legend) {
        p <- p + addLegend(legend.theme)
    } else {
        p <- p
    }

    if (is.null(legend.label$col.point)) {
        p <- p + guides(size = FALSE)
    }

    if (is.null(legend.label$col.branch)) {
        p <- p + guides(color = FALSE)
    }

    p
}

#' Color a branch
#'
#' \code{addBranch} colors a branch or some edges.
#'
#' @param tree A phylo object
#' @param branch A vector of node numbers labels to specify the branches to be
#'   colored. Each branch is represented by its branch node. A leaf node
#'   reprents the edge connecting the leaf and its parent.
#' @param col.branch A vector of colors. Its length should be one or equals to
#'   the length of \strong{branch}. If \strong{col.branch} has the same length
#'   as \strong{branch}, the branches are colored correspondingly with the
#'   \strong{col.branch}. The default is blue.
#' @param col.other A color for the branches other than those specified in
#' \strong{branch}
#' @param addTo NULL or a plot of a phylo object.
#' @param ... see also \code{\link[ggtree]{ggtree}}
#'
#' @import ggplot2
#' @importFrom ggtree ggtree %<+%
#' @return A figure
#' @author Ruizhu Huang
#' @keywords internal
#' @examples
#' # data(tinyTree)
#' # addBranch(tree = tinyTree, branch = 17,
#' # col.branch = "blue", col.other = "grey")

addBranch <- function(tree, branch, col.branch,
                      col.other, addTo = NULL, ...) {

    # node number required
    if (inherits(branch, "character")) {
        branch <- transNode(tree = tree, input = branch,
                            message = FALSE)
    } else {
        branch <- branch
    }

    # -------------------------------------------------------
    # create a data frame to indicate the selected edges
    # -------------------------------------------------------

    p <- ggtree(tree)
    d <- p$data[, "node", drop = FALSE]

    # The edges selected to be colored
    eList <- lapply(branch, findOS, tree = tree,
                    only.Tip = FALSE, self.include = TRUE,
                    return = "number")
    el <- unlist(lapply(eList, length))
    eList <- eList[order(el, decreasing = TRUE)]
    if (length(col.branch) == length(branch)) {
        col.branch <- col.branch[order(el, decreasing = TRUE)]
    }
    dList <- mapply(function(x, y) {
        cbind.data.frame(node = y, group = x,
                         stringsAsFactors = FALSE)},
        x = col.branch, y = eList, SIMPLIFY = FALSE,
        USE.NAMES = FALSE)
    df <- do.call(rbind, dList)

    Truth <- rep("grp_other", nrow(d))
    Truth[match(df$node, d$node)] <- df$group
    d <- cbind.data.frame(d, Truth = Truth, stringsAsFactors = FALSE)

    # return
    if (is.null(addTo)) {
        fig <- ggtree(tree, ...)
    } else {
        fig <- addTo
    }

    fig %<+% d + aes(colour = Truth)

}

#' Add points to the tree plot
#'
#' \code{addPoint} adds points to a plot of phylogenetic tree.
#'
#' @param tree A phylo object
#' @param point A vector of node numbers or node labels to specify the
#' locations to add points in the tree
#' @param col.point A color for the \strong{point}. It has length equal to one.
#' @param addTo NULL or a plot of a phylo object.
#' @param ... see also \code{\link[ggtree]{ggtree}}
#'
#' @import ggplot2
#' @importFrom ggtree ggtree geom_point2
#' @return A figure
#' @author Ruizhu Huang
#' @keywords internal
#' @examples
#' data(tinyTree)
#'
#'
addPoint <- function(tree, point, col.point,
                     addTo = NULL, ...) {
    p <- ggtree(tree)
    d <- p$data[, "node", drop = FALSE]

    # node number required
    if (inherits(point, "character")) {
        point <- transNode(tree = tree, input = point,
                           message = FALSE)
    } else {
        point <- point
    }

    # -------------------------------------------------------
    # create a data frame to store the information for points
    # -------------------------------------------------------
    Estimate <- ifelse(d$node %in% point, "YES_Found",
                       "NO_Found")
    show <- ifelse(d$node %in% point, TRUE, FALSE)
    d <- cbind.data.frame(d, Estimate = Estimate, show = show)

    if (is.null(addTo)) {
        fig <- ggtree(tree, ...)
    } else {
        fig <- addTo
    }

    fig %<+% d +
        geom_point2(aes(subset = show, color = Estimate,
                        size = Estimate))

}

#' Visualize the phylogenetic tree
#'
#' \code{addZoom} zooms in a phylogenetic tree.
#'
#' @param tree A phylo object
#' @param zoomNode A vector of nodes to be zoomed in. If default (NULL), the
#'   tree is not zoomed in.
#' @param zoomLevel A numeric vector. Its length is equal to 1 or equal to the
#'   length of \strong{zoomNode}. If default (NULL), a leaf is zoomed in its
#'   direct parent level and an internal node is zoomed in its own level.
#'
#' @param zoomScale A numeric vector. Its length is equal to one or equal to the
#'   length of \strong{zoomNode}. If \strong{zoomScale} has the same length as
#'   \strong{zoomNode}, the branches are zoomed in with different scales
#'   corresponding to the value of \strong{zoomScale}. If default (NULL), tree
#'   is not zoomed in.
#' @param addTo NULL or a plot of a phylo object.
#' @param ... see also \code{\link[ggtree]{ggtree}}
#'
#' @import ggplot2
#' @importFrom ggtree ggtree %>% scaleClade
#' @return A figure
#' @author Ruizhu Huang
#' @keywords internal
#' @examples
#' # data(tinyTree)
#' # addZoom(tree = tinyTree, zoomNode = 17,
#' # zoomScale = 3)

addZoom <- function(tree, zoomNode = NULL, zoomLevel = NULL,
                    zoomScale = NULL, addTo = NULL, ...) {

    # node number required
    if (is.character(zoomNode)) {
        zoomNode <- transNode(tree = tree, input = zoomNode,
                              message = FALSE)
    } else {
        zoomNode <- zoomNode
    }

    zList <- lapply(zoomNode, findOS, tree = tree,
                    only.Tip = FALSE, self.include = TRUE,
                    return = "number")
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
    nodZ <- findAncestor(tree = tree, node = zoomNode,
                         level = zoomLevel,
                         return = "number")
    names(nodZ) <- names(zoomScale) <- zoomNode
    # remove nodes which are the descendants of the others
    nodZW <- rmDesc(node = nodZ, tree = tree)
    zoomScale[!names(zoomScale) %in% names(nodZW)] <- 1

    if (is.null(addTo)) {
        fig <- ggtree(tree, ...)
    } else {
        fig <- addTo
    }

    # zoom the selected nodes
    i <- 1
    repeat {
        fig <- fig %>% scaleClade(nodZ[i], scale = zoomScale[i])
        i <- i + 1
        if (i > length(nodZ)) {
            break
        }
    }


    lim <- c(min(fig$data$y), max(fig$data$y))

    ## I reset limits for y because ggtree function use ylim to limit y axis.
    ## This would lead to issues, like points not displayed when zoom in some
    ## branches at the case that layout is circular or radical.
    suppressMessages(fig <- fig + scale_y_continuous(limits = lim))

    fig
}

#' Add legend
#' \code{addLegend} customizes the legend.
#'
#' @param legend.theme A list of arguments used for the theme in ggplot2 package
#'   (see \code{\link[ggplot2]{theme}} ) and starting with "legend."
#'
#' @import ggplot2
#' @importFrom utils modifyList
#' @return a list
#' @author Ruizhu Huang
#' @keywords internal


addLegend <- function(legend.theme) {

    # default way to put legend
    li1 <- list(legend.position = "right",
                legend.text = element_text(size= 12),
                legend.key.size = unit(4,"cm"),
                legend.key.height = unit(0.4,"cm"),
                legend.key.width = unit(0.5, "cm"),
                legend.title = element_text(size = 15)
                #,
                # legend.background = element_rect(),
                #legend.box.background = element_rect()
    )
    # user defined
    if (is.null(legend.theme)) {
        legend.theme <- list(NULL)
    }
    li2 <- legend.theme
    # overwrite the default
    li <- modifyList(li1, li2)
    # ggplot2 theme
    do.call(theme, li)
}


#' Customize the scale
#'
#' \code{sizeScale} customizes the size scale.
#'
#' @param col.point A color for the \strong{point}. It has length equal to one.
#' @param size.point The size for the \strong{point}. It has length equal to
#'   one.
#' @param legend.label A list with three members: "col.branch", "col.other", and
#'   "col.point". The elements order in each member matches with the
#'   corresponding argument \strong{col.branch}, \strong{col.other} and
#'   \strong{col.point}, and will display in the legend. See Examples.
#' @param legend.title A vector to specify the title of the legend. It must be
#'   named with "branch" and "point" to match with the argument \strong{branch}
#'   and \strong{point}.
#' @param size.point.legend the point size shown in the legend for
#' \strong{point}.
#' @param legend TRUE or FALSE
#'
#' @import ggplot2
#' @importFrom utils modifyList
#' @return ggproto object (Scale)
#' @author Ruizhu Huang
#' @keywords internal


sizeScale <- function(col.point, size.point,
                      legend.label, legend.title,
                      size.point.legend, legend) {
    # if legend is required, correct the label with guide_legend
    if (legend) {
        ll <- list("branch" = NULL, "point" = NULL)
        lt <- as.list(legend.title)
        names(lt) <- names(legend.title)
        legend.title <- modifyList(ll, lt)
        scale_size_manual(values = size.point,
                          labels = legend.label$col.point,
                          guide = guide_legend(
                              title = legend.title$point,
                              override.aes = list(
                                  shape = 16, color = col.point,
                                  size = size.point.legend)))
    } else {
        scale_size_manual(values = size.point)
    }
}

#' Customize the color
#'
#' \code{colScale} customizes the color scale.
#'
#' @param branch A vector of node numbers labels to specify the branches to be
#'   colored. Each branch is represented by its branch node. A leaf node
#'   reprents the edge connecting the leaf and its parent.
#' @param point A vector of node numbers or node labels to specify the locations
#'   to add points in the tree.
#' @param col.branch A vector of colors. Its length should be one or equals to
#'   the length of \strong{branch}. If \strong{col.branch} has the same length
#'   as \strong{branch}, the branches are colored correspondingly with the
#'   \strong{col.branch}. The default is blue.
#' @param col.other A color for the branches other than those specified in
#'   \strong{branch}
#' @param col.point A color for the \strong{point}. It has length equal to one.
#' @param legend.label A list with three members: "col.branch", "col.other", and
#'   "col.point". The elements order in each member matches with the
#'   corresponding argument \strong{col.branch}, \strong{col.other} and
#'   \strong{col.point}, and will display in the legend. See Examples.
#' @param legend.title A vector to specify the title of the legend. It must be
#'   named with "branch" and "point" to match with the argument \strong{branch}
#'   and \strong{point}.
#' @param size.line.legend the line size shown in the legend for \strong{branch}
#' @param legend TRUE or FALSE. Default is FALSE. If TRUE, the legend is
#'   created.
#'
#' @import ggplot2
#' @importFrom utils modifyList tail
#' @importFrom stats setNames
#' @return ggproto object (color)
#' @author Ruizhu Huang
#' @keywords internal
#'

colScale <- function(branch,
                     point,
                     col.branch,
                     col.other,
                     col.point,
                     legend.label,
                     legend.title,
                     size.line.legend,
                     legend) {
    # colG is to correct the label
    # colV is to correct the value
    # if(is.null(point)){
    #   colG <- colV <- c(col.branch, col.other)
    #   names(colV) <- c(col.branch, "grp_other")
    # }else{
    #   colG <- colV <- c(col.branch, col.other, col.point)
    #   names(colV) <- c(col.branch, "grp_other", "YES_Found")
    # }

    # colG is created to correct the color
    # vG is created to output the label
    if (length(legend.label$col.branch) > length(col.branch)) {
        stop("Same color with different labels. You probably need more colors")
    }

    if (is.null(point)) {
        cG <- list(col.branch, col.other)
        names(cG) <- c("col.branch", "col.other")
        colV <- c(col.branch, col.other)
        names(colV) <- c(col.branch, "grp_other")
    } else {
        cG <- list(col.branch, col.other, col.point)
        names(cG) <- c("col.branch", "col.other", "col.point")
        colV <- c(col.branch, col.other, col.point)
        names(colV) <- c(col.branch, "grp_other", "YES_Found")
    }

    if (legend) {
        #if legend label is not provided
        if (is.null(legend.label)) {
            stop("legend.label isn't provided")
        }

        # decide the content in the legend (branch, other or point)
        # ll is a template
        ll <- list(col.branch = "",
                   col.other = "",
                   col.point = "")
        listG <- listLab <- ll[names(ll) %in% names(cG)]
        listG <- modifyList(listG, cG)
        listLab <- modifyList(listLab, legend.label)

        # check whether listG and listLab have the same composition pattern.
        llG <- lapply(listG, FUN = function(x){match(x, unique(x))})
        llLab <- lapply(listLab, FUN = function(x){match(x, unique(x))})
        if(!setequal(llG, llLab)){
            message("\n The legend label isn't correctly specified. \n")
        }
        # match the color and the label
        namG <- mapply(function(x, y) {
            names(x) <- y
            x
        }, x = setNames(listG, NULL),
        y = setNames(listLab, NULL))

        if (is.list(namG)) {
            colG <- unlist(namG)
        } else {
            colG <- namG
        }

        colG <- colG[!(duplicated(colG) &
                           duplicated(names(colG)))]
        lab <- names(colG)
        ww <- tail(which(lab %in% legend.label$col.point),1)
        # lab <- ifelse(names(colG) %in% legend.label$col.point, "",
        #names(colG))

        lab[ww] <- ""
        #colG <- unlist(setNames(namG, NULL))

        # if there are duplicates for the pairs of color and lable,
        # remove them.
        #colG <- colG[!(duplicated(colG) & duplicated(names(colG)))]
        # lab <- ifelse(names(colG) %in% legend.label$col.point, "",
        #names(colG))
        #lab <- names(colG)
        #lab[3] <- ifelse(lab[3] %in% legend.label$col.point,
        #                "", lab[3] )
        lty <- ifelse(lab %in% "", "blank", "solid")
        du <- duplicated(colG) & duplicated(names(colG))
        lab <- ifelse(du, "", lab)
        lty <- ifelse(du, "blank", lty)


        # update legend.title
        ll <- list("branch" = NULL, "point" = NULL)
        lt <- as.list(legend.title)
        names(lt) <- names(legend.title)
        legend.title <- modifyList(ll, lt)

        scale_color_manual(
            values = colV,
            labels = lab,
            guide = guide_legend(
                title = legend.title$branch,
                override.aes = list(
                    color = colG,
                    linetype = lty,
                    shape = rep(NA, length(colG)),
                    size = size.line.legend
                )
            )
        )
    } else {
        scale_color_manual(values = colV)
    }

}

