#' Sphere surface
#' @description Surface of a sphere (arbitrary dimension).
#'
#' @param d dimension of the space
#' @param r radius of the sphere
#'
#' @return The surface of the sphere of radius \code{r} in the
#' \code{d}-dimensional space.
#' @export
#'
#' @examples
#' r <- 2
#' surface_sphere(3, r)
#' 4*pi*r^2
#' # perimeter of the unit circle:
#' surface_sphere(2)
surface_sphere <- function(d, r=1){
  r^(d-1) * 2 * pi^(d/2)/gamma(d/2)
}

#' Sphere patch surface
#' @description Surface of a sphere patch.
#'
#' @param r radius
#' @param phi1,phi2 numbers defining the latitudinal angle range
#' @param theta1,theta2 numbers defining the longitudinal angle range
#'
#' @return The surface of the sphere patch.
#' @export
#'
#' @seealso \code{\link{surface_stri}} for the surface of a spherical triangle.
#' @details A sphere patch is the part of the sphere whose polar angles
#' \code{theta} and \code{phi} satisfy
#' \code{0 <= theta1 <= theta <= theta2 <= 2*pi} and
#' \code{0 <= phi1 <= phi <= phi2 <= pi}.
#'
#' @examples
#' # surface of the first orthant:
#' surface_spherePatch(r=1, phi1=0, phi2=pi/2, theta1=0, theta2=pi/2)
#' surface_stri(r=1, c(1,0,0), c(0,1,0), c(0,0,1))
surface_spherePatch <- function(r, phi1, phi2, theta1, theta2){
  r*r*(theta2-theta1)*(cos(phi1)-cos(phi2))
}


stri_angles2area <- function(r, angles){
  r * r * ( sum(angles) - pi )
}

stri_sides2angles <- function(r, sides){
  asu <- sides[1L] / r
  bsu <- sides[2L] / r
  csu <- sides[3L] / r
  ssu <- ( asu + bsu + csu ) / 2
  tan_a2 <- sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) /
                     ( sin ( ssu ) * sin ( ssu - asu )     ) )
  tan_b2 <- sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) /
                     ( sin ( ssu ) * sin ( ssu - bsu )     ) )
  tan_c2 <- sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) /
                     ( sin ( ssu ) * sin ( ssu - csu )     ) )
  2 * atan(c(tan_a2, tan_b2, tan_c2))
}

stri_vertices2sides <- function(r, v1, v2, v3){
  r * acos(c(crossprod(v2,v3), crossprod(v3,v1), crossprod(v1,v2)))
}

#' Spherical triangle surface
#' @description Surface of a spherical triangle.
#'
#' @param r radius
#' @param v1,v2,v3 vertices
#'
#' @return The surface of the spherical triangle of radius \code{r} with
#' vertices \code{v1}, \code{v2}, \code{v3}.
#' @export
#'
#' @examples
#' # surface of the first orthant:
#' surface_stri(r=1, c(1,0,0), c(0,1,0), c(0,0,1))
surface_stri <- function(r, v1, v2, v3){
  stri_angles2area(r, stri_sides2angles(r, stri_vertices2sides(r, v1, v2, v3)))
}


#' Triangle surface
#' @description Surface of a triangle.
#'
#' @param v1,v2,v3 vertices of the triangle
#'
#' @return The surface of the triangle with vertices \code{v1}, \code{v2},
#' \code{v3}.
#' @export
#'
#' @examples
#' surface_triangle(c(0,0), c(0,1), c(1,0))
surface_triangle <- function(v1, v2, v3){
  0.5*abs((v2[1L]-v1[1L])*(v3[2L]-v1[2L]) - (v3[1L]-v1[1L])*(v2[2L]-v1[2L]))
}


#' Torus surface
#' @description Surface of a torus.
#' 
#' @param R major radius
#' @param r minor radius
#' 
#' @return The surface area of the torus.
#' @export
surface_torus <- function(R, r){
  4 * pi*pi * R * r
}

#' Spherical cap surface
#' @description Surface of a spherical cap.
#' 
#' @param r radius of the sphere
#' @param h height of the cap
#' 
#' @return The surface area of the spherical cap.
#' @export
surface_sphericalCap <- function(r, h){
  2 * pi * r * h
}