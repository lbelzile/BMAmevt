#' Likelihood function (spectral density) and random generator in the Pairwise Beta and NL models.
#'
#' Applies to angular data sets. The density is given with respect to the Lebesgue measure on \eqn{R^{p-1}}{R^(p-1)}, where \code{p} is the number of columns in \code{x} (or the length of \code{x}, if the latter is a single point).
#' @title Pairwise Beta (PB) and Nested Asymmetric Logistic (NL) distributions
#' @inheritParams prior.pb
#' @param x An angular data set (may be reduced to a single point).
#' A \code{npoints*dimData} matrix
#' (or a vector of length(\code{dimData}).
#' For the NL model, \code{dimData} is always \eqn{3}. 
#' Each row is a point on the simplex, so that the sum of each rows
#' should equal \eqn{1} (the error tolerance is set to \code{1e-8}
#' in this package).
#' @param par The parameter for the Pairwise Beta or the Nested Logistic  density.
#' \itemize{
#' \item In the Pairwise Beta model, \code{par} is of length
#' \code{choose(p,2)+1}. The first element is the global dependence
#' parameter, the subsequent ones are the pairwise dependence
#' parameters, in lexicographic order (\emph{e.g.}
#' \eqn{\beta_{12}, \beta_{13}, \beta_{23}}).
#' \item In the NL model, \code{par} is  a vector of length four  with components between zero and one. The first one is the global dependence parameter, the three subsequent ones are the pairwise dependence parameters, again in lexicographic order. 
#' }
#' @param log Logical. Should the density be returned on the log scale ?
#' @param vectorial Logical.
#' Should a vector or a single value be returned ?
#' @return The value returned by the likelihood function is imposed (see
#' \emph{e.g.} \code{\link{posteriorMCMC}}.
#' In contrast, the random variable have unconstrained output format.
#' \itemize{
#' \item \code{dpairbeta} returns the likelihood  as a  single number if \code{vectorial ==FALSE}, or as a vector of size \code{nrow(x)} containing the likelihood of each angular data point.  If \code{log == TRUE},  the log-likelihood is returned instead.
#' \code{rpairbeta} returns a matrix with \code{n}
#' rows and \code{dimData} columns.
#' \item \code{dnestlog} returns the likelihood  as a  single number if \code{vectorial ==FALSE}, or as a vector of size \code{nrow(x)} containing the likelihood of each angular data point.  If \code{log == TRUE},  the log-likelihood is returned instead.
#' \code{rnestlog} returns a matrix with \code{n} rows and \code{dimData} columns if \code{return.points==FALSE} (the default). Otherwise,
#' a list is returned, with two elements:
#' \itemize{
#'   \item  \code{Angles}: The angular data set
#'   \item   \code{Points}: The full tri-variate data set above
#' \code{threshold} (\emph{i.e.} \code{Angles}
#' multiplied by the radial components)
#'   }
#' }
#' @export
#' @keywords distribution models multivariate
dpairbeta<-
  function(x , 
            par=c(1,rep(2,choose(4,2)+1)), 
           log = FALSE, vectorial = TRUE)
{

  ##evaluates the pairwise beta density for a matrix of  points
  ##(one row = one point)
  ##on the simplex with given pairwise parameters for a
  ##given dimension p , with respect to the lebesgue measure 
#####!!!!!!!!on the projected unit simplex  #####


  xvect = as.double(as.vector(t(x) ))
  
  if(is.vector(x))
    { p = as.integer( length(x) ) 
      n = as.integer(1)
    }  else
    {
      p = as.integer(ncol(x) )
      n = as.integer(nrow(x) )
    }

  if(length(par) <2)
      stop("misspecified parameter")

  alpha= as.double (par[1])
  beta = as.double (par[-1])

  if(vectorial)
      {
        result=as.double(rep(0,n))
      }  else
    {
      result=as.double(0)
    }

  C.out = .C("d_pairbeta", xvect, n , p , alpha, beta, as.integer(log),
    as.integer(vectorial),is_error =as.integer(0),result =  result)
##browser()
  if (C.out$is_error==1)
      {
    ##   browser()
        stop("in d_pairbeta: misbehaved memory allocation")
      }
return(C.out$result)
}
 
