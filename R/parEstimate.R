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
#' @author Ruizhu Huang
#' @examples
#' \dontrun{
#' library(GUniFrac)
#' data("throat.otu.tab")
#' throat.par <- parEstimate(data = t(throat.otu.tab))
#' }

parEstimate <- function(data) {

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
