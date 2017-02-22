\name{transf.to.equi}
\alias{transf.to.equi}
\alias{transf.to.rect}
\title{Linear coordinate transformations}
\usage{
  transf.to.equi(vect)

  transf.to.rect(vect)
}
\arguments{
  \item{vect}{a bi-variate vector, giving the first two
  coordinates of the angular point to be transformed.}
}
\value{
  The vector obtained by linear transformation.
}
\description{
  Switching coordinates system between equilateral and
  right-angled representation of the two-dimensional
  simplex.
}
\details{
  If \code{transf.to.rect}, is called, \code{vect} must
  belongs to the triangle \eqn{[(0,0), (\sqrt(2), 0),
  (\sqrt(2)/2,\sqrt(3/2) ) ]}{[(0,0), (\sqrt(2), 0),
  (\sqrt(2)/2,\sqrt(3/2) ) ]} and the result lies in
  \eqn{([(0,0), (1,0), (0,1)]}{([(0,0), (1,0), (0,1)]}.
  \code{transf.to.equi} is the reciprocal.
}
\examples{
\dontrun{ transf.to.equi(c(sqrt(2)/2, sqrt(3/8) ) )}
}
\author{
  Anne Sabourin
}

