\name{ddirimix.grid}
\alias{ddirimix.grid}
\title{Plots the Dirichlet mixture density on a discretization grid}
\usage{
  ddirimix.grid(par = get("dm.expar.D3k3"), wei = par$wei,
    Mu = par$Mu, lnu = par$lnu, npoints = 30,
    eps = 10^(-3), equi = TRUE, marginal = TRUE,
    coord = c(1, 2, 3), invisible = TRUE, displ = TRUE,
    ...)
}
\arguments{
  \item{marginal}{logical. If \code{TRUE}, the angular
  density corresponds to the marginal intensity measure,
  over coordinates \code{coord}. Otherwise, it is only the
  projection of the full dimensional angular measure (hence
  the moments constraints is not satisfied anymore).}

  \item{coord}{A vector of size 3: the indices of the
  coordinates upon which the marginalization is to be
  done.}

  \item{invisible}{Logical: should the result be returned
  as invisible ?}

  \item{displ}{Logical: should a plot be issued ?}

  \item{...}{Additional arguments to be passed to
  \code{\link{dgridplot}}.}

  \item{par}{The parameter list for the Dirichlet mixture
  model.}

  \item{wei}{Optional. If present, overrides the value of
  \code{par$wei}.}

  \item{Mu}{Optional. If present, overrides the value of
  \code{par$Mu}.}

  \item{lnu}{Optional. If present, overrides the value of
  \code{par$lnu}.}

  \item{npoints}{The number of grid nodes on the squared
  grid containing the desired triangle.}

  \item{eps}{Positive number: minimum distance from any
  node inside the simplex to the simplex boundary}

  \item{equi}{logical. Is the simplex represented as an
  equilateral triangle (if \code{TRUE}) or a right triangle
  (if \code{FALSE}) ?}
}
\value{
  The discretized density
}
\description{
  Only valid in the tri-variate case
}

