\name{dgridplot}
\alias{dgridplot}
\title{Image and/or Contour plots of  spectral densities in  trivariate  extreme value  models}
\usage{
  dgridplot(density = matrix(5 * sin(1/73 * (1:(40 * 40)))^2, ncol = 40, nrow = 40),
    eps = 10^(-3), equi = TRUE, add = FALSE,
    breaks = seq(-0.01, 5.1, length.out = 1000),
    levels = seq(0, 6, length.out = 13),
    col.lines = "black", labcex = 0.8, background = FALSE,
    col.polygon = gray(0.5), lab1 = "w1", lab2 = "w2",
    lab3 = "w3", ...)
}
\arguments{
  \item{density}{A \code{npoints*npoints} matrix containing
  the density's values scattered on the discretization grid
  defined by \code{npoints, equi, eps} (see
  \code{\link{discretize}}).}

  \item{add}{Logical. Should the contours be added to a
  currently active plot ?}

  \item{breaks}{Set of breakpoints for the gray scale
  colors. See \code{\link[graphics]{image}}}

  \item{levels}{Levels to which plot the contour lines. See
  \code{\link[graphics]{contour}}}

  \item{labcex}{\code{cex} for contour labeling. See
  \code{\link[graphics]{contour}}.}

  \item{col.lines}{The color to be used for the contour
  lines.}

  \item{background}{Logical. Should a the background be
  filled inside the simplex \emph{via} a call to
  \code{\link[graphics]{image}} ?}

  \item{...}{Additional graphical parameters and arguments
  to be passed to \code{\link[graphics]{contour}} and
  \code{\link[graphics]{image}}.}

  \item{eps}{Positive number: minimum distance from any
  node inside the simplex to the simplex boundary}

  \item{equi}{logical. Is the simplex represented as an
  equilateral triangle (if \code{TRUE}) or a right triangle
  (if \code{FALSE}) ?}

  \item{col.polygon}{The background color outside the
  simplex.}

  \item{lab1}{Character string: label for first component.}

  \item{lab2}{Character string: label for second
  component.}

  \item{lab3}{Character string: label for third component.}
}
\description{
  Plots contours or gray-scale level sets of a spectral
  density on the two-dimensional simplex.
}
\details{
  The function interprets the \code{density} matrix as
  \code{\link[graphics]{contour}} does, \emph{i.e.} as a
  table of \code{f(X[i], Y[j])} values, with column 1 at
  the bottom, where \code{X} and \code{Y} are returned by
  \code{\link{discretize}} and \code{f} is the density
  function.
}
\examples{
wrapper <- function(x, y, my.fun,...)
      {
       sapply(seq_along(x), FUN = function(i) my.fun(x[i], y[i],...))
      }

grid <- discretize(npoints=40,eps=1e-3,equi=FALSE)

Density <- outer(grid$X,grid$Y,FUN=wrapper,
                my.fun=function(x,y){10*((x/2)^2+y^2)*((x+y)<1)})

dgridplot(density= Density,npoints=40, equi=FALSE)
}

