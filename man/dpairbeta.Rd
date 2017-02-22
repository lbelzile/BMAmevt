\name{dnestlog}
\alias{dnestlog}
\alias{dpairbeta}
\alias{rnestlog}
\alias{rpairbeta}
\title{Pairwise Beta (PB) and Nested Asymmetric Logistic (NL) distributions}
\usage{
  dnestlog(x = rbind(c(0.1, 0.3, 0.6), c(0.3, 0.3, 0.4)),
    par = c(0.5, 0.5, 0.2, 0.3), log = FALSE,
    vectorial = T)

  dpairbeta(x, par = c(1, rep(2, choose(4, 2) + 1)),
    log = FALSE, vectorial = TRUE)

  rnestlog(n = 5, par = c(0.2, 0.3, 0.4, 0.5),
    threshold = 1000, return.points = FALSE)

  rpairbeta(n = 1, dimData = 3, par = c(1, rep(1, 3)))
}
\arguments{
  \item{x}{An angular data set (may be reduced to a single
  point). A \code{npoints*dimData} matrix (or a vector of
  length(\code{dimData}). For the NL model, \code{dimData}
  is always \eqn{3}. Each row is a point on the simplex, so
  that the sum of each rows should equal \eqn{1} (the error
  tolerance is set to \code{1e-8} in this package).}

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

  \item{log}{Logical. Should the density be returned on the
  log scale ?}

  \item{vectorial}{Logical. Should a vector or a single
  value be returned ?}

  \item{threshold}{The radial threshold \eqn{r} above which
  the simulated points should be kept to build the angular
  dataset. Should be set to a high value, for the
  asymptotic approximation \deqn{P(W \in B |\; ||X||
  >r)\simeq H(B)}{P(W \in B | ||X|| >r) ~ H(B)} to hold.}

  \item{return.points}{logical: should the censored
  vectorial dataset corresponding to the angular one be
  returned ?}

  \item{n}{The number of points on the simplex to be
  generated.}

  \item{dimData}{the dimension of the sample space, which
  is \eqn{1 + } the dimension of the simplex.}
}
\value{
  The value returned by the likelihood function is imposed
  (see \emph{e.g.} \code{\link{posteriorMCMC}}. In
  contrast, the random variable have unconstrained output
  format. \itemize{ \item \code{dpairbeta} returns the
  likelihood as a single number if \code{vectorial
  ==FALSE}, or as a vector of size \code{nrow(x)}
  containing the likelihood of each angular data point.  If
  \code{log == TRUE}, the log-likelihood is returned
  instead. \code{rpairbeta} returns a matrix with \code{n}
  rows and \code{dimData} columns. \item \code{dnestlog}
  returns the likelihood as a single number if
  \code{vectorial ==FALSE}, or as a vector of size
  \code{nrow(x)} containing the likelihood of each angular
  data point.  If \code{log == TRUE}, the log-likelihood is
  returned instead. \code{rnestlog} returns a matrix with
  \code{n} rows and \code{dimData} columns if
  \code{return.points==FALSE} (the default). Otherwise, a
  list is returned, with two elements: \itemize{ \item
  \code{Angles}: The angular data set \item \code{Points}:
  The full tri-variate data set above \code{threshold}
  (\emph{i.e.} \code{Angles} multiplied by the radial
  components) } }
}
\description{
  Likelihood function (spectral density) and random
  generator in the Pairwise Beta and NL models.
}
\details{
  Applies to angular data sets. The density is given with
  respect to the Lebesgue measure on
  \eqn{R^{p-1}}{R^(p-1)}, where \code{p} is the number of
  columns in \code{x} (or the length of \code{x}, if the
  latter is a single point).
}
\keyword{datagen}
\keyword{distribution}
\keyword{models}
\keyword{multivariate}

