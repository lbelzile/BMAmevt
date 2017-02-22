\name{rdirichlet}
\alias{rdirichlet}
\title{Dirichlet distribution: random generator}
\usage{
  rdirichlet(n = 1, alpha)
}
\arguments{
  \item{n}{Number of draws}

  \item{alpha}{Dirichlet parameter: a vector of positive
  number}
}
\value{
  A matrix with \code{n} rows and \code{length(alpha)}
  columns
}
\description{
  Dirichlet distribution: random generator
}

