#' Uniform sampling in a tetrahedron
#' @description Uniform sampling in a tetrahedron (in dimension 3).
#'
#' @param n number of simulations
#' @param v1,v2,v3,v4 vertices of the tetrahedron
#'
#' @return The simulations in a \code{n} times \code{3} matrix.
#' @export
#'
#' @seealso \code{\link{runif_in_simplex}} for sampling in a simplex in
#'   arbitrary dimension.
#'
#' @examples
#' \donttest{library(rgl)
#' tetrahedron <- tetrahedron3d()
#' shade3d(tetrahedron, color = "red", alpha = 0.3)
#' vs <- tetrahedron$vb[1L:3L, ]
#' sims <- runif_in_tetrahedron(100, vs[, 1], vs[, 2], vs[, 3], vs[, 4])
#' points3d(sims)}
runif_in_tetrahedron <- function(n, v1, v2, v3, v4){
  out <- matrix(NA_real_, nrow = n, ncol = 3L)
  c <- matrix(runif(n*3), nrow = n, ncol = 3L)
  for(j in 1L:n){
    c1 <- c[j, 1L]
    c2 <- c[j, 2L]
    c3 <- c[j, 3L]
    if( 1 < c1 + c2 ) {
      c1 <- 1 - c1
      c2 <- 1 - c2
    }
    if( 1 < c2 + c3 ) {
      t <- c3
      c3 <- 1 - c1 - c2
      c2 <- 1 - t
    } else if( 1 < (s <- c1 + c2 + c3) ) {
      t <- c3
      c3 <- s - 1
      c1 <- 1 - c2 - t
    }
    c4 <- 1 - (c1 + c2 + c3)

    out[j, ] <- c1 * v1 + c2 * v2 + c3 * v3 + c4 * v4
  }
  out
}
