% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posteriorMean.r
\name{posteriorMean}
\alias{posteriorMean}
\title{Posterior predictive density on the simplex, for three-dimensional extreme value  models.}
\usage{
posteriorMean(
  post.sample,
  FUN = function(par, ...) {
     par
 },
  from = NULL,
  to = NULL,
  thin = 50,
  displ = TRUE,
  ...
)
}
\arguments{
\item{post.sample}{A posterior sample as returned by \code{\link{posteriorMCMC}}}

\item{FUN}{a parameter functional returning a vector.}

\item{from}{Integer or \code{NULL}. If \code{NULL}, the default value is used. Otherwise,  should be greater than \code{post.sample$Nbin}. Indicates the  index where the averaging process should start. Default to \code{post.sample$Nbin +1}}

\item{to}{Integer or \code{NULL}. If \code{NULL}, the default
value is used. Otherwise, must be lower than \code{Nsim+1}.
Indicates  where the averaging process should stop.
Default to \code{post.sample$Nsim}.}

\item{thin}{Thinning interval.}

\item{displ}{logical. Should a plot be produced ?}

\item{...}{Additional parameters to be passed to \code{FUN}.}
}
\value{
A list made of \describe{
\item{values}{A matrix : each column is the result of \code{FUN} applied to a parameter from the posterior sample.}
\item{est.mean}{The posterior mean}
\item{est.sd}{The posterior standard deviation }
}
}
\description{
Computes an approximation of the posterior mean of a parameter functional, based on a  posterior parameters sample.
}
\details{
Only a sub-sample is used: one out of \code{thin} parameters is used
(thinning). Further, only the parameters produced between time
\code{from} and time \code{to} (included) are kept.
}
\seealso{
\code{\link{posteriorMCMC}}.
}
