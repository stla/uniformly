#' Uniform sampling on/in ellipsoid
#' @description Uniform sampling on an ellipsoid or in an ellipsoid. 
#'   The sampling \emph{in} an ellipsoid is available in arbitrary
#'   dimension. The sampling \emph{on} an ellipsoid is available only in 
#'   dimension 2 or 3.
#'
#' @name runif_ellipsoid
#' @rdname runif_ellipsoid
#'
#' @param n number of simulations
#' @param A symmetric positive-definite matrix defining the ellipsoid (see
#'   Details), of size 2 for \code{runif_on_ellipse} and size 2 or 3 for 
#'   \code{runif_on_ellipsoid} (for size 2 these are the same functions)
#' @param r "radius" (see Details)
#'
#' @return The simulations in a matrix with \code{n} rows.
#'
#' @details The ellipsoid is the set of vectors \code{x} satisfying
#'   \code{t(x) \%*\% A \%*\% x == r^2}. For example, for an axis-aligned 
#'   ellipse with horizontal radius \code{a} and vertical radius \code{b}, take 
#'   \code{A=1/diag(c(a^2,b^2))} and \code{r=1}.
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
#' sims <- runif_on_ellipse(60, A, r)
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
runif_on_ellipse <- function(n, A, r){
  stopifnot(r > 0)
  if(!is.matrix(A) || nrow(A) != 2L || ncol(A) != 2L) {
    stop(
      "The matrix `A` should be square of size 2."
    )
  }
  S <- A / (r * r)
  stopifnot(isSymmetric(S))
  e <- eigen(S, symmetric = TRUE)
  eigenvalues <- e$values
  if(any(eigenvalues <= 0)) stop("`A` is not positive.")
  ainv2 <- eigenvalues[1L]; binv2 <- eigenvalues[2L]
  a2 <- 1/ainv2; b2 <- 1/binv2
  a <- sqrt(a2); b <- sqrt(b2)
  if(isTRUE(all.equal(a, b))) {
    return(runif_on_sphere(n, 2L, a))
  }
  rotMatrix <- e$vectors
  ssims <- runif_on_sphere(n, 2L, r = 1)
  allsims <- cbind(a*ssims[, 1L], b*ssims[, 2L])
  p <- a * sqrt(binv2*sq(ssims[, 2L]) + ainv2*sq(ssims[, 1L]))
  accept <- runif(n) < p 
  sims <- allsims[accept, ]
  nsims <- length(accept)
  while(nsims < n) {
    newn <- n - nsims
    ssims <- runif_on_sphere(newn, 2L, r = 1)
    allsims <- cbind(a*ssims[, 1L], b*ssims[, 2L])
    p <- a * sqrt(binv2*sq(ssims[, 2L]) + ainv2*sq(ssims[, 1L]))
    accept <- runif(n) < p 
    sims <- rbind(sims, allsims[accept, ])
    nsims <- nsims + length(accept)
  }
  t(rotMatrix %*% t(sims))
}

#' @rdname runif_ellipsoid
#' @export
runif_on_ellipsoid <- function(n, A, r){
  if(is.matrix(A) && nrow(A) == 2L && ncol(A) == 2L) {
    return(runif_on_ellipse(n, A, r))
  }
  stopifnot(r > 0)
  if(!is.matrix(A) || nrow(A) != 3L || ncol(A) != 3L) {
    stop(
      "This function should be only used for sampling on the boundary of an ", 
      "ellipse or a triaxial ellipsoid only."
    )
  }
  S <- A / (r * r)
  stopifnot(isSymmetric(S))
  e <- eigen(S, symmetric = TRUE)
  eigenvalues <- e$values
  if(any(eigenvalues <= 0)) stop("`A` is not positive.")
  ainv2 <- eigenvalues[1L]; binv2 <- eigenvalues[2L]; cinv2 <- eigenvalues[3L]
  a2 <- 1/ainv2; b2 <- 1/binv2; c2 <- 1/cinv2
  a <- sqrt(a2); b <- sqrt(b2); c <- sqrt(c2)
  if(isTRUE(all.equal(a, c))) {
    return(runif_on_sphere(n, 3L, a))
  }
  rotMatrix <- e$vectors
  ssims <- runif_on_sphere(n, 3L, r = 1)
  allsims <- cbind(a*ssims[, 1L], b*ssims[, 2L], c*ssims[, 3L])
  p <- a * sqrt(
    binv2*sq(ssims[, 2L]) + cinv2*sq(ssims[, 3L]) + ainv2*sq(ssims[, 1L])
  )
  accept <- runif(n) < p 
  sims <- allsims[accept, ]
  nsims <- length(accept)
  while(nsims < n) {
    newn <- n - nsims
    ssims <- runif_on_sphere(newn, 3L, r = 1)
    allsims <- cbind(a*ssims[, 1L], b*ssims[, 2L], c*ssims[, 3L])
    p <- a * sqrt(
      binv2*sq(ssims[, 2L]) + cinv2*sq(ssims[, 3L]) + ainv2*sq(ssims[, 1L])
    )
    accept <- runif(n) < p 
    sims <- rbind(sims, allsims[accept, ])
    nsims <- nsims + length(accept)
  }
  t(rotMatrix %*% t(sims))
}

#' @rdname runif_ellipsoid
#' @export
runif_in_ellipsoid <- function(n, A, r){
  U <- chol(A)
  x <- runif_in_sphere(n, d=ncol(A), r=r)
  t(backsolve(U, t(x)))
}
