#' Uniform sampling in an annulus
#' @description Uniform sampling in an annulus (dimension 2).
#'
#' @param n number of simulations
#' @param O center of the annulus
#' @param r1 inner radius
#' @param r2 outer radius
#'
#' @return The simulations in a \code{n} times \code{2} matrix.
#' @export
#'
#' @examples
#' sims <- runif_in_annulus(100, c(0, 0), 1, 2)
#' plot(sims, xlim = c(-2, 2), ylim = c(-2, 2), asp = 1, pch = 19)
runif_in_annulus <- function(n, O, r1, r2){
  theta <- runif(n, 0, 2*pi)
  v <- runif(n)
  r <- sqrt ( ( 1-v ) * r1*r1 + v * r2*r2 )
  cbind(O[1L] + r*cos(theta), O[2L] + r*sin(theta))
}
