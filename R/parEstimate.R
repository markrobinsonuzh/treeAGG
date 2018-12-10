
# If data is a matrix, the estimated parameters pi and theta is organized as a
# list and returned as output.
# If data is a list, we confirme this list includes two elements pi and theta,
# and directly return data as output
parEstimate.A <- function(data) {

    if (is.list(data)) {
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

#' Parameter estimation in Dirichlet-multinomial distribution
#'
#' \code{parEstimate} is a wrapper function generated from the function
#' \code{\link[dirmult]{dirmult}} with default setting on \code{init},
#' \code{initscalar}, \code{epsilon}, \code{trace} and \code{mode}. It allows
#' the input \code{data} to accept \code{matrix} or
#' \code{leafSummarizedExperiment} and output the value \code{pi} and
#' \code{theta}.
#'
#' @param data A matrix or leafSummarizedExperiment. Samples in the columns and
#'   entities in the rows.
#'
#' @importFrom dirmult dirmult
#' @importFrom methods is
#' @export
#' @return A list including \dQuote{pi} and \dQuote{theta}
#'
#' @author Ruizhu Huang
#' @examples
#'
#' library(S4Vectors)
#' set.seed(1)
#' y <- matrix(rnbinom(100,size=1,mu=10),nrow=10)
#' colnames(y) <- paste("S", 1:10, sep = "")
#' rownames(y) <- tinyTree$tip.label
#'
#'
#' toy_lse <- leafSummarizedExperiment(tree = tinyTree,
#'                                     assays = list(y))
#' res <- parEstimate(data = toy_lse)
#'
#' metadata(res)$assays.par
#'
parEstimate <- function(data) {

    stopifnot(class(data) %in% c("matrix", "list",
                                 "treeSummarizedExperiment"))

    if (is.matrix(data)) {
        out <- parEstimate.A(data = data)
    }

    if (is.list(data)) {
        out <- parEstimate.B(data = data)
    }

    if (is(data, "treeSummarizedExperiment")) {
        out <- parEstimate.C(data = data)
    }
        return(out)
}
