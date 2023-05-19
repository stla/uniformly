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
#' library(uniformly)
#' set.seed(666L)
#' # ellipse parameters
#' A <- rbind(c(2, 1), c(1, 1))
#' r <- 2
#' # plot the ellipse
#' x1 <- seq(-2.5, 2.5, length.out = 100)
#' x2 <- seq(-3, 3, length.out = 100)
#' z <- outer(
#'   x1, x2, FUN = Vectorize(function(x1, x2) t(c(x1, x2)) %*% A %*% c(x1, x2))
#' )
#' contour(x1, x2, z, nlevels = 1, levels = r^2, asp = 1, drawlabels = FALSE)
#' # simulations on the perimeter
#' sims <- runif_on_ellipsoid(60, A, r)
#' points(sims, pch = 19, col = "blue")
#' # simulations in the area
#' sims <- runif_in_ellipsoid(100, A, r)
#' points(sims, pch = 19, col = "green")
#' \donttest{# 3D example ####
#' A <- matrix(c(5,1,1, 1,3,1, 1,1,1), ncol = 3L)
#' r <- 2
#' # draw the ellipsoid
#' library(misc3d)
#' x <- seq(-1, 1, length.out = 50)
#' y <- seq(-1.5, 1.5, length.out = 50)
#' z <- seq(-2.7, 2.7, length.out = 50)
#' g <- as.matrix(expand.grid(x = x, y = y, z = z))
#' voxel <- 
#'   array(apply(g, 1L, function(v) t(v) %*% A %*% v), dim = c(50, 50, 50))
#' isosurface <- computeContour3d(voxel, max(voxel), r^2, x = x, y = y, z = z)
#' drawScene.rgl(makeTriangles(isosurface, alpha = 0.3))
#' # simulate and plot points on ellipsoid
#' library(rgl)
#' sims <- runif_on_ellipsoid(300, A, r)
#' points3d(sims)}
NULL

#' @rdname runif_ellipsoid
#' @export
runif_on_ellipsoid <- function(n, A, r){
  stopifnot(r > 0)
  S <- A / (r * r)
  stopifnot(isSymmetric(S))
  e <- eigen(S, symmetric = TRUE)
  if(any(e$values <= 0)) stop("`A` is not positive.")
  radiisq <- 1 / e$values
  radii <- sqrt(radiisq)
  rot <- e$vectors
  xyz <- t(vapply(radii, function(a) {
    rnorm(n, 0, a)
  }, FUN.VALUE = numeric(n))) 
  d <- sqrt(colSums(xyz*xyz / radiisq))
  sims0 <- t(xyz) / d
  t(rot %*% t(sims0))
}

#' @rdname runif_ellipsoid
#' @export
runif_in_ellipsoid <- function(n, A, r){
  U <- chol(A)
  x <- runif_in_sphere(n, d=ncol(A), r=r)
  t(backsolve(U, t(x)))
}
