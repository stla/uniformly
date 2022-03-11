#' Uniform sampling on/in sphere
#' @description Uniform sampling on a sphere or in a sphere, in arbitrary
#' dimension.
#' @name runif_sphere
#' @rdname runif_sphere
#'
#' @param n number of simulations
#' @param d dimension of the space
#' @param r radius of the sphere
#'
#' @return The simulations in a \code{n} times \code{d} matrix.
#' @importFrom stats runif rnorm
#'
#' @examples
#' sims <- runif_on_sphere(20, d = 2)
#' plot(sims, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, pch = 19)
#' sims <- runif_in_sphere(100, d = 2)
#' plot(sims, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, pch = 19)
NULL

#' @rdname runif_sphere
#' @export
runif_on_sphere <- function(n, d, r=1){
  sims <- matrix(rnorm(n*d), nrow=n, ncol=d)
  r * sims / sqrt(apply(sims, 1L, crossprod))
}

#' @rdname runif_sphere
#' @export
runif_in_sphere <- function(n, d, r=1){
  r * runif_on_sphere(n, d, runif(n)^(1/d))
}


#' Uniform sampling on a spherical patch
#' @description Uniform sampling on a spherical patch (in dimension 3).
#'
#' @param n number of simulations
#' @param r radius
#' @param phi1,phi2 numbers defining the latitudinal angle range
#' @param theta1,theta2 numbers defining the longitudinal angle range
#'
#' @return The simulations in a \code{n} times \code{3} matrix.
#' @export
#'
#' @seealso \code{\link{runif_on_stri}} for sampling on a spherical triangle.
#' @details A sphere patch is the part of the sphere whose polar angles
#' \code{theta} and \code{phi} satisfy
#' \code{0 <= theta1 <= theta <= theta2 <= 2*pi} and
#' \code{0 <= phi1 <= phi <= phi2 <= pi}.
#'
#' @examples
#' # sampling on the first orthant:
#' sims <- 
#'   runif_on_spherePatch(100, phi1 = 0, phi2 = pi/2, theta1 = 0, theta2 = pi/2)
#' \dontrun{
#' library(rgl)
#' spheres3d(0, 0, 0, color = "red", alpha = 0.5)
#' points3d(sims)}
runif_on_spherePatch <- function(n, r=1, phi1, phi2, theta1, theta2){
  sims1 <- runif(n)
  sims2 <- runif(n)
  theta <- theta1*(1-sims1) + theta2*sims1
  cosphi <- cos(phi1)*(1-sims2) + cos(phi2)*sims2
  sinphi <- sin(acos(cosphi))
  unname(r * cbind(cos(theta)*sinphi, sin(theta)*sinphi, cosphi))
}

#' Uniform sampling on a spherical triangle
#' @description Uniform sampling on a spherical triangle (in dimension 3).
#'
#' @param n number of simulations
#' @param r radius
#' @param v1,v2,v3 vertices
#'
#' @return The simulations in a \code{n} times \code{3} matrix.
#' @export
#'
#' @examples
#' # sampling on the first orthant:
#' sims <- runif_on_stri(100, v1 = c(1, 0, 0), v2 = c(0, 1, 0), v3 = c(0, 0, 1))
#' \dontrun{
#' library(rgl)
#' spheres3d(0, 0, 0, color = "red", alpha = 0.5)
#' points3d(sims)}
runif_on_stri <- function(n, r=1, v1, v2, v3){
  sides <- stri_vertices2sides( 1, v1, v2, v3 )
  angles <- stri_sides2angles ( 1, sides )
  area <- stri_angles2area ( 1, angles )
  alpha <- angles[1L]

  out <- matrix(NA_real_, nrow=n, ncol=3L)
  u1 <- runif(n); u2 <- runif(n)

  for(j in 1L:n){
    # Select the new area.
    area_hat = u1[j] * area
    # Compute the sine and cosine of the angle phi.
    s <- sin ( area_hat - alpha )
    t <- cos ( area_hat - alpha )
    # Compute the pair that determines beta_hat.
    u <- t - cos ( alpha )
    v <- s + sin ( alpha ) * cos ( sides[3L] )
    # Q is the cosine of the new edge length b_hat.
    q <- ( ( v * t - u * s ) * cos ( alpha ) - v ) /
      ( ( v * s + u * t ) * sin ( alpha ) )
    # V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
    w <- c(crossprod( v3, v1 ))
    v31 <- v3 - w*v1
    v31 <- v31 / sqrt(c(crossprod(v31)))
    # V4 is the third vertex of the subtriangle V1, V2, V4.
    v4 <- q * v1 + sqrt(1-q*q) * v31
    # Select cos theta, which will sample along the edge from V2 to V4.
    z <- 1.0 - u2[j] * ( 1 - c(crossprod( v4, v2 )) )
    # V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
    w = c(crossprod( v4, v2 ))
    v42 <- v4 - w*v2
    v42 <- v42 / sqrt(c(crossprod(v42)))
    # Construct the point.
    out[j,] <- z * v2 + sqrt(1-z*z)*v42
  }

  r * out
}

#' Uniform sampling on a spherical cap
#' @description Uniform sampling on a spherical cap (in dimension 3).
#'
#' @param n number of simulations
#' @param r radius of the sphere
#' @param h height of the cap
#'
#' @return The simulations in a \code{n} times \code{3} matrix.
#' @export
#'
#' @examples
#' sims <- runif_on_sphericalCap(500, r = 2, h = 1)
#' \dontrun{
#' library(rgl)
#' spheres3d(0, 0, 0, radius = 2, color = "red", alpha = 0.5)
#' points3d(sims)}
runif_on_sphericalCap <- function(n, r = 1, h){
  stopifnot(h > 0, h < 2*r)
  xy <- runif_in_sphere(n, 2L, 1)
  k <- h * apply(xy, 1L, crossprod)
  s <- sqrt(h * (2*r - k))
  cbind(s*xy, r-k)
}

#' Sampling on hemisphere
#' @description Sampling on a hemisphere according to the Phong density
#' (dimension 3).
#'
#' @param n number of simulations
#' @param alpha parameter of the Phong density, a positive number;
#'   \code{0} for uniform sampling (default)
#' @param r radius
#'
#' @return The simulations in a \code{n} times \code{3} matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' library(rgl)
#' sims <- rphong_on_hemisphere(400, alpha = 10)
#' spheres3d(0, 0, 0, color = "red", alpha = 0.5)
#' points3d(sims)}
rphong_on_hemisphere <- function(n, alpha = 0, r = 1){
  cosphi <- runif(n)^(1/(1+alpha))
  sinphi <- sin(acos(cosphi))
  theta <- runif(n, 0, 2*pi)
  r * cbind(cos(theta)*sinphi, sin(theta)*sinphi, cosphi)
}
