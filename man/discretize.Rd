\name{discretize}
\alias{discretize}
\title{Discretization grid builder.}
\usage{
  discretize(npoints = 40, eps = 0.001, equi = FALSE)
}
\arguments{
  \item{npoints}{The number of grid nodes on the squared
  grid containing the desired triangle.}

  \item{eps}{Positive number: minimum distance from any
  node inside the simplex to the simplex boundary}

  \item{equi}{logical. Is the simplex represented as an
  equilateral triangle (if \code{TRUE}) or a right triangle
  (if \code{FALSE}) ?}
}
\value{
  A list containing two elements: \code{X} and \code{Y},
  vectors of size \code{npoints}, the Cartesian coordinates
  of the grid nodes.
}
\description{
  Builds a discretization grid covering the two-dimensional
  unit simplex, with specified number of points and minimal
  distance from the boundary.
}
\details{
  The \code{npoints*npoints} grid covers either the
  equilateral representation of the simplex, or the right
  angled one. In any case, the grid is \emph{rectangular}:
  some nodes lie outside the triangle. Density computations
  on such a grid should handle the case when the point
  passed as argument is outside the simplex (typically, the
  function should return zero in such a case).
}
\note{
  In case \code{equi==TRUE}, \code{epsilon} is the minimum
  distance from any node inside the simplex to the simplex
  boundary, \emph{after transformation} to the right-angled
  representation.
}

