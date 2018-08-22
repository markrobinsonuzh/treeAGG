
# If data is a matrix, the estimated parameters pi and theta is organized as a
# list and returned as output.
# If data is a list, we confirme this list includes two elements pi and theta,
# and directly return data as output
parEstimate.A <- function(data) {

    if (inherits(data, "list")) {
        ind <- setequal(names(data), c("pi", "theta"))
        if (!ind) {
            stop("Error: data is a list;
           it should provide pi and theta")
        }
        parList <- data
    } else {
        DirMultOutput <- dirmult::dirmult(data = t(data))
        # tip proportion
        estP <- DirMultOutput$pi
        names(estP) <- names(DirMultOutput$pi)

        # parameter alpha for dirichlet distribution
        theta <- DirMultOutput$theta
        parList <-  list(pi = estP, theta = theta)
    }

    return(parList)
}

# If data is a list, we confirme this list includes two elements pi and theta,
# and directly return data as output
parEstimate.B <- function(data) {


        ind <- setequal(names(data), c("pi", "theta"))
        if (!ind) {
            stop("Error: data is a list;
                 it should provide pi and theta")
        }
        parList <- data


    return(parList)
}

# If data is a leafSummarizedExperiment, the estimated parameters pi and theta
# is organized as a list and store in metadata with name assays.par
parEstimate.C <- function(data) {

        # node label
        nodeLab <- rowData(data)$nodeLab
        if (is.null(nodeLab)) {
            nodeLab <- rownames(data)
        }

        # estimate parameters
        pars <- metadata(data)$assays.par
        ind <- setequal(names(pars), c("pi", "theta"))

        if (!ind) {
            DirMultOutput <- dirmult::dirmult(data = t(assays(data)[[1]]))
            # tip proportion
            estP <- DirMultOutput$pi
            names(estP) <- nodeLab

            # parameter alpha for dirichlet distribution
            theta <- DirMultOutput$theta
            metadata(data)$assays.par <-  list(pi = estP, theta = theta)
        }

        return(data)

}

#' Estimate parameters for Dirichlet distribution
#'
#' \code{parEstimate} estimates parameters \eqn{\pi} and \eqn{\theta} from
#' a real data (count table). The entity counts within each sample is assumed to
#' follow a multinomial distribution with probabilities,
#' \eqn{ \pi_{1}, \pi_{2}, ..., \pi_{k-1}}. When \eqn{\pi} follows a Dirichlet
#' distribution,  the marginal distribution of X, the entity count, is
#' Dirichlet-multinomial with parameters \eqn{\pi} and \eqn{\theta}.
#' More details see \code{\link[dirmult]{dirmult}}
#'
#' @param data A count table. Samples in the column and entities in the row.
#'
#' @importFrom dirmult dirmult
#' @export
#'
#' @details \code{parEstimate} is created based on
#' \code{\link[dirmult]{dirmult}}).
#'
#' @return A list including \dQuote{pi} and \dQuote{theta}
#' @name parEstimate
#' @author Ruizhu Huang
#' @examples
#' \dontrun{
#' library(GUniFrac)
#' data("throat.otu.tab")
#' throat.par <- parEstimate(data = t(throat.otu.tab))
#' }
#'
setGeneric("parEstimate", function(data) {
    standardGeneric("parEstimate")
})

#' @rdname parEstimate
setMethod("parEstimate", signature(data = "matrix"), parEstimate.A)

#' @rdname parEstimate
setMethod("parEstimate", signature(data = "list"), parEstimate.B)


#' @rdname parEstimate
setMethod("parEstimate", signature(data = "leafSummarizedExperiment"),
          parEstimate.C)
