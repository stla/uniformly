#' Uniform sampling on/in a triangle
#' @description Uniform sampling on or in a triangle (dimension 2).
#'
#' @name runif_triangle
#' @rdname runif_triangle
#'
#' @param n number of simulations
#' @param v1,v2,v3 vertices of the triangle
#'
#' @return The simulations in a \code{n} times \code{2} matrix.
#'
#' @examples
#' sims <- runif_on_triangle(30, c(0,0), c(1,0), c(0,1))
#' plot(sims, xlim = c(0,1), ylim = c(0,1), pch = 19)
#' sims <- runif_in_triangle(100, c(0,0), c(1,0), c(0,1))
#' plot(sims, xlim = c(0,1), ylim = c(0,1), pch = 19)
NULL

#' @rdname runif_triangle
#' @export
runif_in_triangle <- function(n, v1, v2, v3){
  r1 <- runif(n)
  r2 <- sqrt(runif(n))
  a <- 1 - r2
  b <- r2 * (1-r1)
  c <- r1 * r2
  cbind(a*v1[1L]+b*v2[1L]+c*v3[1L], a*v1[2L]+b*v2[2L]+c*v3[2L])
}

#' @rdname runif_triangle
#' @export
runif_on_triangle <- function(n, v1, v2, v3){
  l1 <- sqrt(c(crossprod(v2-v1)))
  l2 <- sqrt(c(crossprod(v3-v2)))
  l3 <- sqrt(c(crossprod(v1-v3)))
  r <- runif(n, 0, l1+l2+l3)
  out <- matrix(NA_real_, nrow=n, ncol=2L)
  for(j in 1L:n){
    if( r[j] <= l1 ){ # Case 1: between V1 and V2
      s <- ( l1 - r[j] ) / l1
      t <-        r[j]   / l1
      out[j,] <- s*v1 + t*v2
    }else if( r[j] <= l1 + l2 ){ # Case 2: between V2 and V3
      s = ( l2 - r[j] + l1 ) / l2
      t = (      r[j] - l1 ) / l2
      out[j,] <- s*v2 + t*v3
    }else{ # Case 3: between V3 and V1
      s = ( l3 - r[j] + l1 + l2 ) / l3
      t = (      r[j] - l1 - l2 ) / l3
      out[j,] <- s*v3 + t*v1
    }
  }
  out
}
