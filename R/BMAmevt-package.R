##'Toolkit  for Bayesian estimation of the dependence structure
##' in Multivariate Extreme Value parametric models, with possible use of Bayesian model Averaging techniques
##'
##' \tabular{ll}{
##' Package: \tab BMAmevt \cr
##' Type: \tab Package\cr
##' Version: \tab 1.0\cr
##' Date: \tab 2012-03-1\cr
##' License: \tab GPL-2\cr
##' LazyData: \tab no\cr
##' } Includes a Generic MC MC sampler. Estimation of the marginal
##' distributions is a prerequisite, \emph{e.g.} using one of the
##' packages
##' \code{ismev}, \code{evd}, \code{evdbayes} or \code{POT}. This package handles data sets which are assumed
##' to be  marginally unit-Frechet distributed. 
##' @name BMAmevt-package
##' @aliases BMAmevt
##' @docType package
##' @title Bayesian Model Averaging for  Multivariate Extremes
##' @author Anne Sabourin
##' \email{sabourin@@math.univ-lyon1.fr}
##' @keywords package models multivariate htest
##' @seealso \code{evdbayes} 
NULL

##' Five-dimensional air quality dataset recorded in Leeds(U.K.), during five winter seasons.
##' 
##' Contains 590  daily  maxima of five air pollutants
##' (respectively PM10, N0, NO2, 03, S02) recorded in  Leeds (U.K.)
##' during five winter seasons  (1994-1998, November-February included).  Contains NA's.
##' @name winterdat
##' @docType data
##' @format A \eqn{590*5} matrix.
##' @source   \url{http://www.airquality.co.uk}
##' @keywords datasets
NULL

##' Multivariate data set with margins following unit Frechet
##' distribution.
##'
##'  Five-variate dataset which margins follow unit-Frechet distributions,
##' obtained from \code{\link{winterdat}} by probability integral
##' transform.
##' Marginal estimation  was performed by maximum likelihood estimation of a Generalized Pareto distribution over marginal thresholds corresponding to \eqn{0.7} quantiles, following 
##' Cooley \emph{et.al.} (see reference below). The \dQuote{non extreme} part of the marginal distributions was approximated by the empirical distribution function. 
##' @name frechetdat
##' @docType data
##' @format A \eqn{601*5} - matrix: 
##' @references COOLEY, D., DAVIS , R. and  NAVEAU , P. (2010). The pairwise beta distribution: A flexible parametric multivariate  model for extremes. \emph{Journal of Multivariate Analysis 101, 2103-2117}.
##' @keywords datasets
NULL
##' Tri-variate \sQuote{angular} data set approximately distributed according to a multivariate extremes angular distribution
##'
##' The data set is constructed from coordinates (columns)  \eqn{1,2,3} of \code{\link{frechetdat}}.
##' It  contains 100 angular points corresponding to the tri-variate vectors \eqn{V=(X,Y,Z)}  with largest \eqn{L^1}{L1} norm (\eqn{||V||=X+Y+Z}). The angular points are obtained by \sQuote{normalizing}: \emph{e.g.},
##' \eqn{x=X/||V||}. Thus, 
##' each row in \code{Leeds} is a point on the two-dimensional simplex :  \eqn{x+y+z=1}.
##' @name Leeds
##' @docType data
##' @format A \eqn{100*3} - matrix. 
##' @references COOLEY, D., DAVIS , R. and  NAVEAU , P. (2010). The pairwise beta distribution: A flexible parametric multivariate  model for extremes. \emph{Journal of Multivariate Analysis 101, 2103-2117}
##' 
##' RESNICK , S. (1987). Extreme values, regular variation, and point processes, \emph{Applied Probability. A, vol. 4, 
##' Series of the Applied Probability Trust. Springer-Verlag, New York}.
##' @keywords datasets
NULL
##' Default hyper-parameters for the Pairwise Beta model. 
##'
##' The log-transformed dependence parameters are a priori independent, Gausian. This list contains the means and standard deviation for the prior distributions.
##' @name pb.Hpar
##' @docType data
##' @format A list of four parameters: \describe{
##' 
##' \item{ mean.alpha}{
##' Mean of the  log-transformed global dependence parameter. Default to \eqn{0}  )
##' }
##'
##' \item{sd.alpha}{Standard deviation  of the log-transformed global dependence parameter.  Default to \eqn{3}.
##' }
##'\item{mean.beta}{
##' Mean of each of the  log-transformed pairwise  dependence parameters.
##' Default to \eqn{0}  )
##' }
##' \item{sd.beta}{Standard deviation  of each of the  log-transformed
##' pairwise dependence parameters.  Default to \eqn{3}.
##' }
##' }
##' @keywords dataset
NULL
##' Default MC MC tuning parameter for the  Pairwise Beta model. 
##'
##' The proposal for the log-transformed parameters are Gaussian, centered at the current value. 
##' @name pb.MCpar
##' @docType data
##' @format A list made of a  single element: \code{sd},
##' the standard deviation of the normal proposition kernel (on the log-transformed parameter). Default to 
##' \eqn{0.35}.
##' @keywords dataset
NULL
##' Default hyper-parameters for the NL model. 
##'
##' The logit-transformed parameters for the NL model are \emph{a priori}
##' Gaussian. The list has the same format as \code{\link{pb.Hpar}}.
##' @name nl.Hpar
##' @docType data
##' @format A list of four parameters: \describe{
##' \item{mean.alpha, sd.alpha}{%
##' Mean and standard deviation of the normal prior distribution for the  logit-transformed global dependence parameter \eqn{alpha} .
##' Default to \eqn{0, 3}.
##' }
##' \item{mean.beta, sd.beta}{Idem for the pairwise dependence parameters.
##' }
##' }
##' @keywords dataset
NULL
##' Default MC MC tuning parameter for the Nested Asymmetric logistic model. 
##'
##' The proposals (on the logit-scale) are Gaussian, centered aroud the current value. 
##' @name nl.MCpar
##' @docType data
##' @format A list made of a  single element: \code{sd}. The standard deviation of the normal proposition kernel centered at the (logit-transformed)
##' current state. Default to \eqn{0.35}.
##' @keywords dataset
NULL

##' Example of valid Dirichlet mixture parameter for tri-variate extremes.
##' 
##' The Dirichlet mixture density has three components, the center of mass  of the three columns of \code{Mu}, with weights \code{wei} is \eqn{(1/3,1/3,1/3)}: the centroid of the two dimensional unit simplex.
##' @name dm.expar.D3k3
##' @docType data
##' @format A list made of \describe{
##' \item{Mu}{ A \eqn{3*3} matrix, which rows sum to one, such that the
##' center of mass of the three column vectors (weighted with \code{wei}) is the centroid of the simplex: each column is the center of a Dirichlet mixture component. }
##' \item{wei}{A vector of length three, summing to one: the mixture weights}
##' \item{lnu}{ A vector of length three: the logarithm of the concentration parameters. }
##' }
##' @keywords dataset
NULL

##' @importFrom grDevices dev.new gray
##' 
NULL



