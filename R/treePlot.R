#' visualize the phylogenetic tree
#'
#' \code{treePlot} is to visualize the phylogenetic tree. It provides arguments to color and zoom in specified branches.
#'
#' @param tree a phylo object
#' @param trueLoc a vector of nodes (node numbers or labels)
#' @param level a numeric vector specified the number of steps one could take from trueLoc to get closer to the root
#' @param zoomScale a scale to zoom in
#' @param col.trueLoc colors for the branches with root nodes specified in trueLoc. It could be a vector with length one or a vector having the same lenght as trueLoc.
#' @param col.other a color for branches other than those specified by trueLoc
#'
#'
#' @export
#'
#' @return a tree plot
#'
#' @examples
#'
#' data(exTree)
#'
#' treePlot(tree = exTree, trueLoc = c(77, 95, 16),
#' level = c(1, 1, 2), zoomScale = 8,
#' col.trueLoc = c('red','blue','orange'),
#' branch.length = 'none')+
#' geom_text2(aes(label = node))

#'
treePlot <- function(tree, trueLoc, level, zoomScale, col.trueLoc, col.other = "grey",
    col.aggLoc, aggLoc, size.point = 2, size.point.legend = 4, size.line.legend = 2,
    ...) {


    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    # add points to selected nodes
    if (inherits(aggLoc, "character")) {
        aggLoc <- tx_node(tree = tree, input = aggLoc)
    } else {
        aggLoc <- aggLoc
    }

    if (length(col.trueLoc) != length(trueLoc) & length(col.trueLoc) != 1) {
        stop("The length of col.trueLoc should be one or equal to the length of trueLoc")
    }

    # organize the color used in the figure
    if (is.null(names(col.trueLoc))) {
        if (length(col.trueLoc) != length(trueLoc)) {
            names(col.trueLoc) <- "Diff"
        } else {
            names(col.trueLoc) <- trueLoc
        }

    } else {
        col.trueLoc <- col.trueLoc
    }

    if (is.null(names(col.other))) {
        cat("A named vector is requried by col.other; The name is generated automatically as 'Non-diff'")
        names(col.other) <- "Non-diff"
    } else {
        col.other <- col.other
    }

    if (is.null(names(col.aggLoc))) {
        cat("A named vector is requried by col.aggLoc; The name is generated automatically as 'Diff'")
        names(col.aggLoc) <- "Diff"
    } else {
        col.aggLoc <- col.aggLoc
    }


    colG <- c(col.trueLoc, col.other, col.aggLoc)


    p <- ggtree(tree)

    # find the edges selected by the tree aggregation
    desList <- lapply(trueLoc, FUN = function(x) {
        findOS(ancestor = x, tree = tree,
               only.Tip = FALSE, self.include = TRUE)
    })
    names(desList) <- trueLoc
    des_len <- unlist(lapply(desList, length))
    des.1 <- desList[order(des_len, decreasing = TRUE)]
    group_var <- lapply(seq_along(des.1), FUN = function(x) {
                    rep(names(des.1)[x], length(des.1[[x]]))})
    df.1 <- data.frame(node = unlist(des.1, use.names = FALSE),
                       group = unlist(group_var),
                       stringsAsFactors = FALSE)

    if (length(col.trueLoc) == 1) {
        df.1$group <- names(col.trueLoc)
    }

    # create a data frame with three columns: node, Truth, Estimate node (node number);
    # Truth (NO = non-differential, )
    d <- p$data[, "node", drop = FALSE]
    d$Truth <- rep(names(col.other), nrow(d))
    d$Truth[match(df.1$node, d$node)] <- df.1$group
    d$Estimate <- ifelse(d$node %in% aggLoc, names(col.aggLoc), "NO_Found")
    d$show <- ifelse(d$node %in% aggLoc, TRUE, FALSE)



    p1 <- ggtree(tree, aes(colour = Truth), layout = "circular") %<+%
      d +
      geom_point2(aes(subset = show, color = Estimate,
                          size = Estimate)) +
      scale_color_manual(values = colG,
                         labels = ifelse(names(colG) != names(col.aggLoc),
                                         names(colG), ""),
                         guide = guide_legend(
                           title = "Truth",
                           override.aes =list(
                             color = colG, linetype = ifelse(
                               names(colG) == names(col.aggLoc),
                               "blank", "solid"),
                             shape = rep(NA, length(col.trueLoc) + 2),
                             size = size.line.legend))) +
      scale_size_manual(values = size.point,
                        guide = guide_legend(
                          title = "Estimate",
                          override.aes = list(
                            shape = 16, color = col.aggLoc,
                            size = size.point.legend)))


    # define level
    if (is.null(level)) {
        level <- ifelse(des_len > 1, 0, 1)
    } else {
        if (length(level) == 1) {
            level <- rep(level, length(trueLoc))
        } else {
            level <- level
        }
    }
    names(level) <- trueLoc

    # define zoomScale
    if (is.null(zoomScale)) {
        zoomScale <- rep(1, length(trueLoc))
    } else {
        zoomScale <- rep(zoomScale, length(trueLoc))
    }

    # the nodes to be zoomed in
    nodZ <- findAncestor(tree = tree, node = trueLoc, level = level)
    names(nodZ) <- names(zoomScale) <- trueLoc
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

    return(p1)
}
