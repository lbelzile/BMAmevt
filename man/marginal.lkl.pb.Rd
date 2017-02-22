\name{marginal.lkl.nl}
\alias{marginal.lkl.nl}
\alias{marginal.lkl.pb}
\title{Marginal likelihoods of the PB and NL models.}
\usage{
  marginal.lkl.nl(dat, Nsim = 10000, displ = TRUE,
    Hpar = get("nl.Hpar"), Nsim.min = Nsim, precision = 0,
    show.progress = floor(seq(1, Nsim, length.out = 20)))

  marginal.lkl.pb(dat, Nsim = 10000, displ = TRUE,
    Hpar = get("pb.Hpar"), Nsim.min = Nsim, precision = 0,
    show.progress = floor(seq(1, Nsim, length.out = 20)))
}
\arguments{
  \item{dat}{The angular data set relative to which the
  marginal model likelihood is to be computed}

  \item{displ}{logical. If \code{TRUE}, a plot is produced,
  showing the temporal evolution of the cumulative mean,
  with approximate confidence intervals of \eqn{+/-2}
  estimated standard errors.}

  \item{precision}{the desired relative precision. See
  \code{\link{MCpriorIntFun}}.}

  \item{Nsim}{Maximum number of iterations}

  \item{Nsim.min}{The minimum number of iterations to be
  performed.}

  \item{show.progress}{same as in
  \code{\link{posteriorMCMC}}}

  \item{Hpar}{A list containing Hyper-parameters to be
  passed to \code{prior}.}
}
\value{
  The list returned by \code{\link{marginal.lkl}},
  \emph{i.e.}, the one returned by
  \code{\link{MCpriorIntFun}}
}
\description{
  Wrappers for \code{\link{marginal.lkl}}, in the specific
  cases of the PB and NL models, with parameter
  \code{likelihood} set to \code{dpairbeta} or
  \code{dnestlog}, and \code{prior} set to \code{prior.pb}
  or \code{prior.nl}. See \code{\link{MCpriorIntFun}} for
  more details.
}
\examples{
\dontrun{

marginal.lkl.pb(dat=Leeds ,
         Nsim=20e+3 ,
         displ=TRUE, Hpar = get("pb.Hpar") ,
          )

marginal.lkl.nl(dat=Leeds ,
         Nsim=10e+3 ,
         displ=TRUE, Hpar = get("nl.Hpar") ,
          )
}
}
\seealso{
  \code{\link{marginal.lkl}}, \code{\link{MCpriorIntFun}} .
}

