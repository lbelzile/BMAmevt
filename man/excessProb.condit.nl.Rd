\name{excessProb.condit.nl}
\alias{excessProb.condit.nl}
\title{Probability of joint threshold excess in the NL model}
\usage{
  excessProb.condit.nl(par = c(0.3, 0.4, 0.5, 0.6),
    thres = rep(100, 3))
}
\arguments{
  \item{par}{The Nested logistic parameter: of length
  four.}

  \item{thres}{a positive vector of size three.}
}
\value{
  The approximate probability of joint excess, valid when
  at least one coordinate of \code{thres} is large
}
\description{
  Probability of joint threshold excess in the NL model
}

