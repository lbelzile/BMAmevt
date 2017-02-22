\name{dnestlog.grid}
\alias{dnestlog.grid}
\alias{dpairbeta.grid}
\title{PB and NL spectral  densities on the two-dimensional simplex}
\usage{
  dnestlog.grid(par, npoints = 50, eps = 0.001,
    equi = TRUE, displ = TRUE, invisible = TRUE, ...)

  dpairbeta.grid(par, npoints = 50, eps = 0.001,
    equi = TRUE, displ = TRUE, invisible = TRUE, ...)
}
\arguments{
  \item{displ}{logical. Should a plot be produced ?}

  \item{invisible}{logical. If \code{TRUE}, the result is
  returned as \code{invisible}.}

  \item{...}{Additional arguments to be passed to
  \code{\link{dgridplot}}}

  \item{par}{The parameter for the Pairwise Beta or the
  Nested Logistic density. \itemize{ \item In the Pairwise
  Beta model, \code{par} is of length \code{choose(p,2)+1}.
  The first element is the global dependence parameter, the
  subsequent ones are the pairwise dependence parameters,
  in lexicographic order (\emph{e.g.} \eqn{\beta_{12},
  \beta_{13}, \beta_{23}}). \item In the NL model,
  \code{par} is a vector of length four with components
  between zero and one. The first one is the global
  dependence parameter, the three subsequent ones are the
  pairwise dependence parameters, again in lexicographic
  order. }}

  \item{npoints}{The number of grid nodes on the squared
  grid containing the desired triangle.}

  \item{eps}{Positive number: minimum distance from any
  node inside the simplex to the simplex boundary}

  \item{equi}{logical. Is the simplex represented as an
  equilateral triangle (if \code{TRUE}) or a right triangle
  (if \code{FALSE}) ?}
}
\value{
  A \code{npoints*npoints} matrix containing the considered
  density's values on the grid. The row (resp. column)
  indices increase with the first (resp. second) coordinate
  on the simplex.
}
\description{
  The two functions compute respectively the NL and PB
  spectral densities, in the three-dimensional case, on a
  discretization grid. A plot is issued (optional).
}
\note{
  If \code{equi==TRUE}, the density is relative to the
  Hausdorff measure on the simplex itself: the values
  obtained with \code{equi = FALSE} are thus divided by
  \eqn{\sqrt 3}.
}
\examples{
dpairbeta.grid(par=c( 0.8, 8, 5, 2),
npoints=70, eps = 1e-3, equi = TRUE, displ = TRUE, invisible=TRUE)

##  or ...

Dens <- dpairbeta.grid(par=c(0.8, 8, 5, 2),
npoints=70, eps = 1e-3, equi = TRUE, displ = FALSE)
Grid=discretize(npoints=70,eps=1e-3,equi=TRUE)
dev.new()
image(Grid$X, Grid$Y, Dens)
contour(Grid$X, Grid$Y, Dens, add=TRUE)
add.frame(equi=TRUE, npoints=70, axes=FALSE)
}

