#' Uniform sampling in a p-ball
#' @description Uniform sampling in a p-ball (arbitrary dimension).
#'
#' @param n number of simulations
#' @param d dimension
#' @param p exponent in the p-norm, a positive number
#' @param r positive number, the radius
#'
#' @return The simulations in a \code{n} times \code{d} matrix.
#' @export
#' @importFrom pgnorm rpgnorm
#'
#' @examples
#' sims <- runif_in_pball(500, d = 2, p = 1)
#' plot(sims, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
runif_in_pball <- function(n, d, p, r = 1){
  out <- matrix(NA_real_, nrow=n, ncol=d)
  epsilon <- matrix(pgnorm::rpgnorm(n*d, p), nrow=n, ncol=d)
  signs <- matrix(2L*sample.int(2L, n*d, replace=TRUE)-1L, nrow=n, ncol=d)
  x <- signs * epsilon
  z <- runif(n)^(1/d)
  for(j in 1L:n){
    out[j, ] <- r*z[j]*x[j, ]/sum(abs(x[j, ])^p)^(1/p)
  }
  out
}
