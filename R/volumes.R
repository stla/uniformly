#' Sphere volume
#' @description Volume of a sphere (arbitrary dimension).
#'
#' @param d dimension of the space
#' @param r radius of the sphere
#'
#' @return The volume of the sphere with radius \code{r} in the
#' \code{d}-dimensional space.
#' @export
#'
#' @examples
#' r <- 2
#' volume_sphere(3, r)
#' 4/3*pi*r^3
volume_sphere <- function(d, r=1){
  r^d * pi^(d/2)/gamma(d/2 + 1)
}

#' Ellipsoid volume
#' @description Volume of an ellipsoid (arbitrary dimension).
#'
#' @param A symmetric positive-definite matrix defining the ellipsoid (see
#' Details)
#' @param r "radius" (see Details)
#'
#' @return The volume of the ellipsoid.
#' @export
#'
#' @details The (boundary of the) ellipsoid is the set of vectors \code{x}
#' satisfying \code{t(x) \%*\% A \%*\% x == r^2}.
#'
#' @examples
#' # dimension 2 (area), with diagonal matrix A
#' A <- diag(c(2,3))
#' r <- 2
#' volume_ellipsoid(A, r)
#' pi * r^2 / sqrt(A[1,1]*A[2,2])
volume_ellipsoid <- function(A, r){
  volume_sphere(nrow(A),r) / sqrt(det(A))
}

#' Unit simplex volume
#' @description Volume of the unit simplex (arbitrary dimension).
#'
#' @param d dimension of the space
#'
#' @return The volume of the unit simplex in the space of dimension \code{d}.
#'
#' @seealso \code{\link{volume_simplex}} for the volume of an arbitrary simplex.
#'
#' @export
volume_unitSimplex <- function(d){
  1/factorial(d)
}

#' Simplex volume
#' @description Volume of a simplex (arbitrary dimension).
#'
#' @param simplex a \code{(d+1)} times \code{d} matrix giving the vertices of
#' the simplex (rows)
#'
#' @return The volume of the simplex.
#' @export
#'
#' @examples
#' set.seed(666)
#' simplex <- matrix(rnorm(4*3), nrow=4, ncol=3)
#' volume_simplex(simplex)
#' volume_tetrahedron(simplex[1,], simplex[2,], simplex[3,], simplex[4,])
volume_simplex <- function(simplex){
  n <- ncol(simplex)
  stopifnot(is.matrix(simplex), nrow(simplex) == n+1L)
  V <- simplex[-1L,] - matrix(simplex[1L,], nrow = n, ncol = n, byrow = TRUE)
  return(abs(det(V))/factorial(n))
}

extprod3d <- function(x,y){
  c(x[2L]*y[3L] - x[3L]*y[2L],
    x[3L]*y[1L] - x[1L]*y[3L],
    x[1L]*y[2L] - x[2L]*y[1L])
}

#' Tetrahedron volume
#' @description Volume of a tetrahedron (dimension 3).
#'
#' @param v1,v2,v3,v4 vertices of the tetrahedron
#'
#' @return The volume of the tetrahedron.
#' @export
#'
#' @seealso \code{\link{volume_simplex}} for the volume of a simplex in
#' arbitrary dimension.
#'
#' @examples
#' v1 <- c(0,0,0); v2 <- c(1,0,0); v3 <- c(0,1,0); v4 <- c(0,0,1)
#' volume_tetrahedron(v1, v2, v3, v4)
#' volume_unitSimplex(3)
volume_tetrahedron <- function(v1,v2,v3,v4){
  abs(c(crossprod(v1-v4, extprod3d(v2-v4,v3-v4))))/6
}


#' p-ball volume
#' @description Euclidean volume of a p-ball (arbitrary dimension).
#'
#' @param d dimension
#' @param p exponent in the p-norm, a positive number
#' @param r radius of the ball
#'
#' @return The volume of the \code{p}-ball with radius \code{r}.
#' @export
#'
#' @examples
#' volume_pball(d=4, p=2, r=2)
#' volume_sphere(d=4, r=2)
volume_pball <- function(d, p, r=1){
  (2*gamma(1+1/p)*r)^d / gamma(d/p+1)
}


#' Torus volume
#' @description Volume of a torus.
#' 
#' @param R major radius
#' @param r minor radius
#' 
#' @return The volume of the torus.
#' @export
volume_torus <- function(R, r){
  2 * pi*pi * R * r*r
}

#' Spherical cap volume
#' @description Volume of a spherical cap.
#' 
#' @param r radius of the sphere
#' @param h height of the cap
#' 
#' @return The volume of the spherical cap.
#' @export
volume_sphericalCap <- function(r, h){
  pi * h * h * (r - h/3)
}