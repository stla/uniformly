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
  out <- matrix(NA_real_, nrow=n, ncol=3L)
  c <- cbind(matrix(runif(n*3), nrow=n, ncol=3L), NA_real_)
  for(j in 1L:n){
    if( 1 < c[j,1L] + c[j,2L] ){
      c[j,1L] <- 1 - c[j,1L]
      c[j,2L] <- 1 - c[j,2L]
    }
    if( 1 < c[j,2L] + c[j,3L] ){
      t <- c[j,3L]
      c[j,3L] <- 1 - c[j,1L] - c[j,2L]
      c[j,2L] <- 1 - t
    }
    else if( 1 < (s <- sum(c[j,1L:3L])) ){
      t <- c[j,3L]
      c[j,3L] <- s - 1
      c[j,1L] <- 1 - c[j,2L] - t
    }
    c[j,4L] = 1 - sum(c[j,1L:3L])

    out[j,] <- c[j,1L] * v1 + c[j,2L] * v2 + c[j,3L] * v3 + c[j,4L] * v4
  }
  out
}
