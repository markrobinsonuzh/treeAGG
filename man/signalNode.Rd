% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signalNode.R
\name{signalNode}
\alias{signalNode}
\title{find the optimal nodes to short result.}
\usage{
signalNode(tree, node, use.alias = FALSE)
}
\arguments{
\item{tree}{A tree (phylo object)}

\item{node}{A vector of node numbers or node labels}

\item{use.alias}{A logical value, TRUE or FALSE. The default is FALSE, and
the node label would be used to name the output; otherwise, the alias of
node label would be used to name the output. The alias of node label is
created by adding a prefix \code{"Node_"} to the node number if the node is
an internal node or adding a prefix \code{"Leaf_"} if the node is a leaf
node.}
}
\value{
A vector of nodes. The numeric value is the node number, and the
  vector name is the corresponding node label. If a node has no label, it
  would have NA as name when \code{use.alias = FALSE}, and have the alias of
  node label as name when \code{use.alias = TRUE}.
}
\description{
\code{signalNode} is to represent some nodes with their ancestor to make
result as short as possible. The ancestors share exactly the same leaves as
the original nodes.
}
\examples{

data(tinyTree)
library(ggtree)

# PLOT tree
# The node labels are in orange texts and the node numbers are in blue
ggtree(tinyTree,branch.length = 'none')+
    geom_text2(aes(label = label), color = "darkorange",
           hjust = -0.1, vjust = -0.7) +
    geom_text2(aes(label = node), color = "darkblue",
               hjust = -0.5, vjust = 0.7)

## find the node shared by provided node labels
signalNode(node = c('t4','t9'), tree = tinyTree)
signalNode(node = c('t4','t9'), tree = tinyTree)
signalNode(node = c('t10','Node_18', 't8'), tree = tinyTree,
 use.alias = FALSE)
signalNode(node = c('t10','Node_18', 't8'), tree = tinyTree,
 use.alias = TRUE)

## find the node shared by provided node numbers
signalNode(node = c(2, 3), tree = tinyTree)
signalNode(node = c(2, 3, 16), tree = tinyTree)

}
\author{
Ruizhu Huang
}
