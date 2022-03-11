#' Uniform sampling on/in torus
#' @description Uniform sampling on or in a torus (dimension 3).
#'
#' @name runif_torus
#' @rdname runif_torus
#'
#' @param n number of simulations
#' @param R major radius
#' @param r minor radius 
#'
#' @return The simulations in a \code{n} times \code{3} matrix.
#'
#' @examples
#' R <- 3; r <- 2
#' sims_on <- runif_on_torus(50, R = R, r = r)
#' sims_in <- runif_in_torus(50, R = R, r = r)
#' \donttest{library(misc3d)
#' fx <- function(u,v) (R+r*cos(u)) * cos(v)
#' fy <- function(u,v) (R+r*cos(u)) * sin(v)
#' fz <- function(u,v) r*sin(u)
#' parametric3d(
#'   fx, fy, fz, umin = 0, umax = 2*pi, vmin = 0, vmax = 2*pi, alpha = 0.3
#' )
#' library(rgl)
#' points3d(sims_on)
#' points3d(sims_in, color = "red")}
NULL

#' @rdname runif_torus
#' @export
runif_on_torus <- function(n, R, r){
  out <- matrix(NA_real_, nrow=n, ncol=3L)
  j <- 0L
  Phi <- runif(n, 0, 2*pi)
  while(j < n){
    Theta <- runif(1L, 0, 2*pi)
    h <- R + r*cos(Theta)
    if(runif(1L, 0, R+r) <= h){
      j <- j+1L
      out[j,] <- c(h*cos(Phi[j]),  h*sin(Phi[j]), r*sin(Theta))
    }else{
      next
    }
  }
  out
}

#' @rdname runif_torus
#' @export
runif_in_torus <- function(n, R, r){
  out <- matrix(NA_real_, nrow=n, ncol=3L)
  j <- 0L
  while(j < n){
    U <- sqrt(runif(1L))
    alpha <- runif(1L, 0, 2*pi)
    X <- R + r*U*cos(alpha)
    if(runif(1L, 0, R+r) > X) next
    j <- j + 1L
    Theta <- runif(1L, 0, 2*pi)
    out[j, ] <- c(cos(Theta)*X, sin(Theta)*X, r*U*sin(alpha))
  }
  out  
}
