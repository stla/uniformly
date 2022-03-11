#' Uniform sampling in a polygon
#' @description Uniform sampling in a polygon (dimension 2).
#'
#' @param n number of simulations
#' @param vertices two-columns matrix giving the vertices (rows); the vertices
#'   must be ordered (clockwise or counterclockwise)
#' @param center a point with respect to which the polygon is star-shaped, or
#'   \code{"centroid"} (default) to take the centroid (see Details)
#'
#' @return The simulations in a \code{n} times \code{2} matrix.
#' @export
#'
#' @details This function works for a star-shaped polygon, that is, a polygon
#'   that contains a point from which the entire polygon boundary is visible.
#'   This point must be given in the \code{center} argument. If the polygon is
#'   convex, any point inside the polygon is suitable (thus the default option
#'   of the \code{center} argument is appropriate in this case).
#'
#' @examples
#' vs <- matrix(c(0.951056516295154, 0.309016994374947,
#'                0.224513988289793, 0.309016994374947,
#'                -0.951056516295154, 0.309016994374948,
#'                -0.363271264002681, -0.118033988749895,
#'                0.587785252292473, -0.809016994374948,
#'                0.36327126400268, -0.118033988749895,
#'                0, 1,
#'                -0.224513988289793, 0.309016994374947,
#'                -0.587785252292473, -0.809016994374947,
#'                0, -0.381966011250105),
#'              ncol=2, byrow=TRUE)
#' sims <- runif_in_polygon(500, vs)
#' plot(sims, xlim = c(-1, 1), ylim = c(-1, 1), pch = 19, asp = 1)
runif_in_polygon <- function(n, vertices, center = "centroid"){
  out <- matrix(NA_real_, nrow=n, ncol=2L)
  if(identical(center, "centroid")){
    center <- colMeans(vertices)
  }
  nv <- nrow(vertices)
  areas <- numeric(nv)
  for(i in 1L:nv){
    ip1 <- ifelse(i<nv, i+1L, 1L)
    areas[i] <- surface_triangle(vertices[i, ], vertices[ip1, ], center)
  }
  areas <- areas / sum(areas)
  for(j in 1L:n){
    t1 <- sample.int(nv, 1L, prob = areas)
    t2 <- ifelse(t1<nv, t1+1L, 1L)
    out[j, ] <- runif_in_triangle(1L, vertices[t1, ], vertices[t2, ], center)
  }
  out
}
