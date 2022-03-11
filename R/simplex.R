#' Uniform sampling on/in a unit simplex
#' @description Uniform sampling on or in a unit simplex (arbitrary dimension).
#'
#' @name runif_unitSimplex
#' @rdname runif_unitSimplex
#'
#' @param n number of simulations
#' @param d dimension of the space
#'
#' @return The simulations in a \code{n} times \code{d} matrix.
#'
#' @importFrom stats rexp
#'
#' @seealso \code{\link{runif_in_tetrahedron}} for sampling in an arbitrary
#'   tetrahedron in dimension 3; \code{\link{runif_in_simplex}} for sampling
#'   in an arbitrary simplex.
#'
#' @examples
#' \donttest{library(rgl)
#' sims <- runif_on_unitSimplex(300, d = 3)
#' points3d(sims)}
NULL

#' @rdname runif_unitSimplex
#' @export
runif_on_unitSimplex <- function(n, d){
  e <- matrix(rexp(n*d,1), nrow=n, ncol=d)
  out <- e / apply(e, 1L, sum)
  u <- runif(n)
  h <- 1/(1+sqrt(d))
  for(j in 1L:n){
    if(h < u[j]){
      out[j, sample.int(d,1L)] <- 0
    }
  }
  out
}

#' @rdname runif_unitSimplex
#' @export
runif_in_unitSimplex <- function(n, d){
  e <- matrix(rexp(n*(d+1), 1), nrow=n, ncol=d+1L)
  e[, -1L] / rowSums(e)
}


#' Uniform sampling in a simplex
#' @description Uniform sampling in a simplex (arbitrary dimension).
#'
#' @param n number of simulations
#' @param simplex a \code{(d+1)} times \code{d} matrix giving the vertices of
#'   the simplex (rows)
#'
#' @return The simulations in a \code{n} times \code{d} matrix.
#' @export
#'
#' @note In dimension 3, you can use \code{\link{runif_in_tetrahedron}} instead.
#'
#' @examples
#' simplex <- rbind(c(0,0,0), c(1,0,0), c(1,1,0), c(1,1,2))
#' sims <- runif_in_simplex(1000, simplex)
#' \donttest{library(rgl)
#' points3d(sims)}
runif_in_simplex <- function(n, simplex){
  d <- ncol(simplex); nv <- nrow(simplex)
  stopifnot(is.matrix(simplex), nv == d + 1L)
  out <- matrix(NA_real_, nrow=n, ncol=d)
  U <- matrix(runif(n*d), nrow=n, ncol=d)
  for(j in 1L:n){
    p <- simplex[1L,]
    for(i in 2L:nv){
      u <- U[j,i-1L]^(1/(i-1))
      p <- u*p + (1-u)*simplex[i,]
    }
    out[j,] <- p
  }
  out
}
