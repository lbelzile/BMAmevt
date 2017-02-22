\name{excessProb.nl}
\alias{excessProb.nl}
\title{Posterior distribution the probability of joint threshold excess, in the NL model.}
\usage{
  excessProb.nl(post.sample, from = NULL, to = NULL,
    thin = 100, thres = rep(100, 3), known.par = FALSE,
    true.par, displ = FALSE)
}
\arguments{
  \item{post.sample}{The posterior sample, as returned by
  \code{posteriorMCMC}}

  \item{known.par}{logical. Is the true parameter known ?}

  \item{true.par}{The true parameter, only used if
  \code{known.par=TRUE}}

  \item{from}{Integer or \code{NULL}. If \code{NULL}, the
  default value is used. Otherwise, should be greater than
  \code{post.sample$Nbin}. Indicates the index where the
  averaging process should start. Default to
  \code{post.sample$Nbin +1}}

  \item{to}{Integer or \code{NULL}. If \code{NULL}, the
  default value is used. Otherwise, must be lower than
  \code{Nsim+1}. Indicates where the averaging process
  should stop. Default to \code{post.sample$Nsim}.}

  \item{thin}{Thinning interval.}

  \item{displ}{logical. Should a plot be produced ?}

  \item{thres}{a positive vector of size three.}
}
\value{
  A list made of \describe{ \item{whole}{The output of
  \code{posteriorMean} called with
  \code{FUN=excessProb.condit.nl}.} \item{mean}{The
  posterior mean of the excess probability}
  \item{esterr}{The standard deviation of the mean
  estimator} \item{estsd}{The standard deviation of the
  excess probability, in the posterior sample. }
  \item{lowquant}{The lower 0.1 quantile of the empirical
  posterior distribution of the excess probability }
  \item{upquant}{The upper 0.1 quantile of the empirical
  posterior distribution of the excess probability }
  \item{true}{\code{NULL} if \code{known.par=FALSE},
  otherwise the excess probability in the true model.} }
}
\description{
  Posterior distribution the probability of joint threshold
  excess, in the NL model.
}

