#' Uniform sampling on/in ellipsoid
#' @description Uniform sampling on an ellipsoid or in an ellipsoid (arbitrary
#'   dimension).
#'
#' @name runif_ellipsoid
#' @rdname runif_ellipsoid
#'
#' @param n number of simulations
#' @param A symmetric positive-definite matrix defining the ellipsoid (see
#'   Details)
#' @param r "radius" (see Details)
#'
#' @return The simulations in a matrix with \code{n} rows.
#'
#' @details The ellipsoid is the set of vectors \code{x} satisfying
#' \code{t(x) \%*\% A \%*\% x == r^2}.
#'
#' @examples
#' A <- rbind(c(2, 1), c(1, 1))
#' r <- 2
#' sims <- runif_on_ellipsoid(30, A, r)
#' plot(sims, xlim = c(-2, 2), ylim = c(-3, 3), asp = 1, pch = 19)
#' sims <- runif_in_ellipsoid(100, A, r)
#' plot(sims, xlim = c(-2, 2), ylim = c(-3, 3), asp = 1, pch = 19)
#' \donttest{# 3D example
#' A <- matrix(c(5,1,1, 1,3,1, 1,1,1), ncol = 3L)
#' r <- 2
#' # draw the ellipsoid
#' library(misc3d)
#' x <- seq(-1, 1, len = 50)
#' y <- seq(-1.5, 1.5, len = 50)
#' z <- seq(-2.7, 2.7, len = 50)
#' g <- as.matrix(expand.grid(x = x, y = y, z = z))
#' voxel <- 
#'   array(apply(g, 1L, function(v) t(v) %*% A %*% v), dim = c(50, 50, 50))
#' isosurface <- computeContour3d(voxel, max(voxel), r^2, x = x, y = y, z = z)
#' drawScene.rgl(makeTriangles(isosurface, alpha = 0.3))
#' # simulate and plot points on ellipsoid
#' library(rgl)
#' sims <- runif_on_ellipsoid(200, A, r)
#' points3d(sims)}
NULL

#' @rdname runif_ellipsoid
#' @export
runif_on_ellipsoid <- function(n, A, r){
  U <- chol(A)
  x <- runif_on_sphere(n, d=ncol(A), r=r)
  t(backsolve(U, t(x)))
}

#' @rdname runif_ellipsoid
#' @export
runif_in_ellipsoid <- function(n, A, r){
  U <- chol(A)
  x <- runif_in_sphere(n, d=ncol(A), r=r)
  t(backsolve(U, t(x)))
}
