\name{prior.nl}
\alias{prior.nl}
\title{Prior parameter distribution for the NL  model}
\usage{
  prior.nl(type = c("r", "d"), n, par, Hpar, log,
    dimData = 3)
}
\arguments{
  \item{type}{One of the character strings \code{"r"},
  \code{"d"}}

  \item{n}{The number of parameters to be generated. Only
  used if \code{type == "r"}.}

  \item{par}{A vector of length four, with component
  comprised between \eqn{0} and \eqn{1} (both end points
  excluded for the first element and \eqn{1} included for
  the others): The parameter where the density is to be
  taken. Only used if \code{type=="d"}.

  In the NL model, \code{par} is of length \eqn{4}. The
  first element is the global dependence parameter, the
  others are partial dependence parameter between pairs
  (12), (13), (23) respectively.

  In the NL model, \code{par} is of length \eqn{4}. The
  first element has the same interpretation as in the NL
  model, the subsequent ones are dependence parameters
  between}

  \item{Hpar}{list of Hyper-parameters : see
  \code{\link{nl.Hpar}} for a template.}

  \item{log}{logical. Should the density be returned on the
  log scale ? Only used if \code{type=="d"}}

  \item{dimData}{The dimension of the sample space, equal
  to \eqn{3}. Only for compatibility with \emph{e.g.}
  \code{\link{posteriorMCMC}}.}
}
\value{
  Either a matrix with \code{n} rows containing a random
  parameter sample generated under the prior (if type ==
  "d"), or the (log)-density of the parameter \code{par}.
}
\description{
  Density and generating function of the prior
  distribution.
}
\details{
  The four parameters are independent, the
  logit-transformed parameters follow a normal
  distribution.
}
\examples{
\dontrun{prior.nl(type="r", n=5 ,Hpar=get("nl.Hpar")) }

\dontrun{prior.trinl(type="r", n=5 ,Hpar=get("nl.Hpar")) }
\dontrun{prior.pb(type="d", par=rep(0.5,2), Hpar=get("nl.Hpar")) }
}
\author{
  Anne Sabourin
}

