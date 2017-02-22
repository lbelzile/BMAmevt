\name{expfunction.nl}
\alias{expfunction.nl}
\title{Exponent function in the NL  model.}
\usage{
  expfunction.nl(par = c(0.3, 0.4, 0.5, 0.6),
    x = 10 * rep(1, 3))
}
\arguments{
  \item{par}{The parameter for the NL distribution,
  respectively of length two or four.}

  \item{x}{A vector of three extended positive real
  numbers}
}
\value{
  the value of \eqn{V(x)} for \eqn{x=thres}.
}
\description{
  The exponent function \eqn{V} for a max-stable variable
  \eqn{M} is such that \eqn{P(M<x) = exp(-V(x))}
}

