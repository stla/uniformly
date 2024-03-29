% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sphere.R
\name{rphong_on_hemisphere}
\alias{rphong_on_hemisphere}
\title{Sampling on hemisphere}
\usage{
rphong_on_hemisphere(n, alpha = 0, r = 1)
}
\arguments{
\item{n}{number of simulations}

\item{alpha}{parameter of the Phong density, a positive number;
\code{0} for uniform sampling (default)}

\item{r}{radius}
}
\value{
The simulations in a \code{n} times \code{3} matrix.
}
\description{
Sampling on a hemisphere according to the Phong density
(dimension 3).
}
\examples{
\dontrun{
library(rgl)
sims <- rphong_on_hemisphere(400, alpha = 10)
spheres3d(0, 0, 0, color = "red", alpha = 0.5)
points3d(sims)}
}
