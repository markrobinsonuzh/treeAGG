% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viewSim.R
\name{viewSim}
\alias{viewSim}
\title{visualize simulated scenario}
\usage{
viewSim(obj, layout, zoomScale, legend.theme = list(legend.position = c(0.15,
  0.6)), tip.label = FALSE)
}
\arguments{
\item{obj}{The output from \code{simData}}

\item{layout}{The tree layout. Select one from 'rectangular', 'slanted', 'fan', 'circular', 'radial', 'equal_angle' and 'daylight'.}

\item{zoomScale}{A positive numeric value. If it is above one, branches with fold change equal to one (non-signal branch) will be zoomed in; If below one, they
will be shrinked.}

\item{legend.theme}{A list of arguments used for the theme in ggplot2
package (see \code{\link[ggplot2]{theme}} ) and starting with "legend."}

\item{tip.label}{TRUE or FALSE. Default is FALSE. If TRUE, the leaves with
fold change above or below 1 will be labelled.}
}
\value{
a figure
}
\description{
\code{viewSim} is to visualize the output from the function \code{simData}.
}