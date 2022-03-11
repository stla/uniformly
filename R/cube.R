#' Uniform sampling on/in cube
#' @description Uniform sampling on or in a cube (arbitrary dimension).
#'
#' @name runif_cube
#' @rdname runif_cube
#'
#' @param n number of simulations
#' @param d dimension
#' @param O center of the cube
#' @param r radius (half-side) of the cube
#'
#' @return The simulations in a \code{n} times \code{d} matrix.
#'
#' @examples
#' sims <- runif_on_cube(60, d = 2)
#' plot(sims, xlim = c(-1,1), ylim = c(-1,1), pch = 19, asp = 1)
#' \donttest{sims <- runif_in_cube(50, d = 3)
#' library(scatterplot3d)
#' scatterplot3d(sims, pch = 19, highlight.3d = TRUE, asp = 1)}
NULL

#' @rdname runif_cube
#' @export
runif_in_cube <- function(n, d, O = rep(0, d), r = 1){
  sweep(matrix(runif(n*d, -r, r), nrow=n, ncol=d), 2L, O, "+")
}

#' @rdname runif_cube
#' @export
runif_on_cube <- function(n, d, O=rep(0,d), r=1){
  out <- matrix(runif(n*d, -r, r), nrow=n, ncol=d)
  i <- sample.int(d, n, replace = TRUE)
  b <- ifelse(sample.int(2L, n, replace = TRUE)==1L, -r, r)
  for(j in 1L:n){
    out[j, i[j]] <- b[j]
  }
  sweep(out, 2L, O, "+")
}
