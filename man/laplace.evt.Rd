\name{laplace.evt}
\alias{laplace.evt}
\title{Laplace approximation of  a model marginal likelihood by Laplace approximation.}
\usage{
  laplace.evt(mode = NULL, npar = 4, likelihood, prior,
    Hpar, data, link, unlink, method = "L-BFGS-B")
}
\arguments{
  \item{mode}{The parameter vector (on the
  \dQuote{unlinked} scale, \emph{i.e.} before
  transformation to the real line) which maximizes the
  posterior density, or \code{NULL}.}

  \item{npar}{The size of the parameter vector. Default to
  four.}

  \item{likelihood}{The likelihood function, \emph{e.g.}
  \code{\link{dpairbeta}} or \code{\link{dnestlog}}}

  \item{prior}{The prior density (takes an
  \dQuote{unlinked} parameter as argument and returns the
  density of the \code{linked} parameter)}

  \item{Hpar}{The prior hyper parameter list.}

  \item{data}{The angular dataset}

  \item{link}{The link function, from the
  \dQuote{classical} or \dQuote{unlinked} parametrization
  onto the real line. (\emph{e.g.} \code{log} for the PB
  model, an \code{logit} for the NL model)}

  \item{unlink}{The inverse link function (\emph{e.g.}
  \code{exp} for the PB model and \code{invlogit} for the
  NL model)}

  \item{method}{The optimization method to be used. Default
  to \code{"L-BFGS-B"}.}
}
\value{
  A list made of \describe{ \item{mode}{the parameter (on
  the unlinked scale) deemed to maximize the posterior
  density. This is equal to the argument if the latter is
  not null.} \item{value}{The value of the posterior,
  evaluated at \code{mode}.} \item{laplace.llh}{The
  logarithm of the estimated marginal likelihood}
  \item{invHess}{The inverse of the estimated hessian
  matrix at \code{mode}} }
}
\description{
  Approximation of a model marginal likelihood by Laplace
  method.
}
\details{
  The posterior mode is either supplied, or approximated by
  numerical optimization. For an introduction about
  Laplace's method, see \emph{e.g.} Kass and Raftery, 1995
  and the references therein.
}
\references{
  KASS, R.E. and RAFTERY, A.E. (1995). Bayes Factors.
  \emph{Journal of the American Statistical Association,
  Vol. 90, No.430}
}

