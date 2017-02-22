\name{posteriorMCMC}
\alias{posteriorMCMC}
\title{MC MC sampler for parametric spectral measures}
\usage{
  posteriorMCMC(prior = function(type = c("r", "d"), n, par, Hpar, log, dimData) {
        NULL },
    proposal = function(type = c("r", "d"), cur.par, prop.par, MCpar, log) {
        NULL },
    likelihood = function(x, par, log, vectorial) {
        NULL }, Nsim, dat, Hpar, MCpar, Nbin = 0,
    par.start = NULL,
    show.progress = floor(seq(1, Nsim, length.out = 20)),
    seed = NULL, kind = "Mersenne-Twister", save = FALSE,
    class = NULL, name.save = NULL, save.directory = "~",
    name.dat = "", name.model = "")
}
\arguments{
  \item{Nsim}{Total number of iterations to perform.}

  \item{Nbin}{Length of the burn-in period.}

  \item{par.start}{Starting point for the MC MC sampler.}

  \item{dat}{An angular data set, \emph{e.g.} constructed
  by \code{\link{cons.angular.dat}}: A matrix which rows
  are the Cartesian coordinates of points on the unit
  simplex (summing to one).}

  \item{likelihood}{The likelihood function. Should be of
  type\cr \code{function(x, par, log, vectorial)}, where
  \code{log} and \code{vectorial} are logical flags
  indicating respectively if the result is to be returned
  on the log-scale, and if the value is a vector of length
  \code{nrow(x)} or a single number (the likelihood, or the
  log-likelihood, for the data set \code{x}). See
  \code{\link{dpairbeta}} or \code{\link{dnestlog}} for
  templates.}

  \item{proposal}{The proposal function: of type \cr
  \code{function(type = c("r","d"), cur.par, prop.par,
  MCpar, log) }. Should return the (logarithm of) the
  proposal density for the move \code{cur.par --> prop.par}
  if \code{type == "d"}. If \code{type =="r"},
  \code{proposal} must return a candidate parameter,
  depending on \code{cur.par}, as a vector. See
  \code{\link{proposal.pb}} or \code{\link{proposal.nl}}
  for templates.}

  \item{prior}{The prior distribution: of type \cr
  \code{function(type=c("r","d"), n ,par, Hpar, log,
  dimData )}, where \code{dimData} is the dimension of the
  sample space (\emph{e.g.}, for the two-dimensional
  simplex (triangle), \code{dimData=3}. Should return
  either a matrix with \code{n} rows containing a random
  parameter sample generated under the prior (if \code{type
  == "d"}), or the density of the parameter \code{par} (the
  logarithm of the density if \code{log==TRUE}. See
  \code{\link{prior.pb}} and \code{\link{prior.nl}} for
  templates.}

  \item{Hpar}{A list containing Hyper-parameters to be
  passed to \code{prior}.}

  \item{MCpar}{A list containing MC MC tuning parameters to
  be passed to \code{proposal}.}

  \item{show.progress}{An vector of integers containing the
  times (iteration numbers) at which a message showing
  progression will be printed on the standard output.}

  \item{seed}{The seed to be set \emph{via}
  \code{\link[base]{set.seed}}.}

  \item{kind}{The kind of random numbers generator. Default
  to "Mersenne-Twister". See \code{\link[base]{set.seed}}
  for details.}

  \item{save}{Logical. Should the result be saved ?}

  \item{class}{Optional character string: additional class
  attribute to be assigned to the result. A predefined
  class \code{"PBNLpostsample"} exists, for which a method
  performing convergence diagnostics is defined (see
  \code{\link{diagnose}} )}

  \item{save.directory}{A character string giving the
  directory where the result is to be saved (without
  trailing slash).}

  \item{name.save}{A character string giving the name under
  which the result is to be saved. If \code{NULL}
  (default), nothing is saved. Otherwise, the result is
  saved in file \cr \code{paste(save.directory,"/",
  name.save,".rda",sep="")}. A "log" list is also saved,
  named \cr \code{paste(name.save, ".log", sep="")}, in
  file \cr \code{paste(save.directory,"/",
  name.log,".rda",sep="")}.}

  \item{name.dat}{A character string naming the data set
  used for inference. Default to \code{""}.}

  \item{name.model}{A character string naming the model.
  Default to \code{""}.}
}
\value{
  A list made of \itemize{ \item \code{stored.vals}: A
  \code{(Nsim-Nbin)*d} matrix, where \code{d} is the
  dimension of the parameter space. \item \code{llh} A
  vector of size \code{(Nsim-Nbin)} containing the
  loglikelihoods evaluated at each parameter of the
  posterior sample. \item \code{lprior} A vector of size
  \code{(Nsim-Nbin)} containing the logarithm of the prior
  densities evaluated at each parameter of the posterior
  sample. \item \code{elapsed}: The time elapsed, as given
  by \code{proc.time} between the start and the end of the
  run. \item \code{Nsim}: The same as the passed argument
  \item \code{Nbin}: idem. \item\code{n.accept}: The total
  number of accepted proposals. \item \code{n.accept.kept}:
  The number of accepted proposals after the burn-in
  period. \item \code{emp.mean} The estimated posterior
  parameters mean \item \code{emp.sd} The empirical
  posterior sample standard deviation.}
}
\description{
  Generates a posterior parameters sample, and computes the
  posterior mean and component-wise variance on-line.
}
\examples{
data(Leeds)
data(pb.Hpar)
data(pb.MCpar)
postsample1 <- posteriorMCMC(Nsim=1e+3,Nbin=500,
         dat= Leeds,
         prior = prior.pb,
         proposal = proposal.pb,
         likelihood = dpairbeta,
         Hpar=pb.Hpar,
         MCpar=pb.MCpar)

dim(postsample1[[1]])
postsample1[-1]

\dontrun{
## a more realistic one:

postsample2 <- posteriorMCMC(Nsim=50e+3,Nbin=15e+3,
         dat= Leeds,
         prior = prior.pb,
         proposal = proposal.pb,
         likelihood = dpairbeta,
         Hpar=pb.Hpar,
         MCpar=pb.MCpar)
dim(postsample2[[1]])
postsample2[-1]
}
}
\seealso{
  \code{\link{posteriorMCMC.pb}},
  \code{\link{posteriorMCMC.pb}} for specific uses in the
  PB and the NL models.
}
\keyword{htest}
\keyword{multivariate}

