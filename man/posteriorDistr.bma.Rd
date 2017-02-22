\name{posteriorDistr.bma}
\alias{posteriorDistr.bma}
\title{Posterior distribution in the average model}
\usage{
  posteriorDistr.bma(postweights = c(0.5, 0.5),
    post.distrs = list())
}
\arguments{
  \item{postweights}{a vector of positive real numbers,
  summing to one: the posterior weights (in the same order
  as the elements of \code{post.distrs}) of the individual
  models.}

  \item{post.distrs}{A list of same length as
  \code{postweights}. Each element must be a vector which
  will be used as a posterior sample.}
}
\value{
  A matrix with two rows and as many columns as the sum of
  the lengths of the elements of \code{post.distrs}. The
  second line contains the weighted posterior sample in the
  BMA; the first line contains the weights to be assigned
  to each corresponding value on the second one.
}
\description{
  Builds an empirical distribution defined as a sum of
  weighted Dirac masses from posterior samples in
  individual models.
}

