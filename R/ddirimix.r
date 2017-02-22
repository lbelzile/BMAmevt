##' Likelihood function (spectral density on the simplex)
##' and angular data  sampler   in the  Dirichlet mixture model.
##'
##'  The spectral probability measure  defined  on the simplex
##' characterizes the
##' dependence structure of multivariate extreme value models.
##' The parameter list for a mixture
##' with  \eqn{k}  components, is made of
##' \describe{
##' \item{Mu}{ The density kernel centers
##' \eqn{\mu_{i,m}, 1\le i \le p, 1\le m \le k}{\mu[1:p,1:k]} :
##' A  \eqn{p*k} matrix,
##' which columns sum to one, and such that \code{Mu \%*\% wei=1},
##' for the moments constraint to be satisfied. 
##' Each column is a Dirichlet kernel center.
##' }
##' \item{wei}{ The weights vector for the kernel densities:
##' A vector of  \eqn{k} positive numbers summing to one.}
##' \item{lnu}{The logarithms of the shape parameters
##' \eqn{nu_m, 1\le m \le k}{\nu[1:k] } for the density kernels:
##' a vector of size \eqn{k}.}
##' }
##' The moments constraint imposes  that the barycenter of the columns in
##' \code{Mu}, with weights \code{wei}, be the center of the simplex.
##' @title Angular density/likelihood function in the Dirichlet Mixture
##' model.
##' @param x An angular data set which may be reduced to a single point: 
##' A \eqn{n*p} matrix or a vector of length \code{p}, where
##' \eqn{p} is the dimension of the sample space and \eqn{n} is
##' the sample size.
##' Each row is a point on the simplex, so that  each row sum to one. 
##' The error tolerance is set to \code{1e-8}
##' in this package.
##' @param par The parameter list for the Dirichlet mixture model.
##'  @param wei Optional. If present, overrides the value of
##' \code{par$wei}.
##' @param Mu Optional. If present, overrides the value of
##' \code{par$Mu}.
##' @param lnu Optional. If present, overrides the value of
##' \code{par$lnu}.
##' @param log Logical: should the density or the likelihood be returned on the log-scale ?
##' @param vectorial Logical: Should a vector of size \eqn{n} or a single value be returned ? 
##' @return \code{ddirimix} returns the likelihood  as a  single number if
##' \code{vectorial ==FALSE}, or as a vector of size
##' \code{nrow(x)} containing the likelihood of each angular data point.
##' If \code{log == TRUE},  the log-likelihood is returned instead.
##' \code{rdirimix} returns a matrix with \code{n} points and
##' \code{p=nrow(Mu)} columns. 
##' @export
ddirimix <-
function(x=c(0.1,0.2,0.7),
           par, #= get("dm.expar.D3k3"),
           wei=par$wei,
           Mu= par$Mu,
           lnu=par$lnu,
           log=FALSE,
           vectorial=FALSE)
  {
    k=length(lnu)
    nu=exp(lnu)
    if(is.vector(x))
      {
        x=matrix(x,nrow=1)
      }
    n=nrow(x)
    p=ncol(x)

    if(vectorial)
      {result=as.double(rep(0,n))}
    else
      {
        result=as.double(0)
      }
    
    C.out=.C("d_dirimix",as.double(t(x)),as.integer(n),as.integer(k),
      as.integer(p),as.double(wei),as.double(Mu),as.double(nu),
      as.integer(log),as.integer(vectorial),is_error=as.integer(0),
      result=result)

    if (C.out$is_error==1)
      {
        stop("in ddirimix: some data or kernels not on the simplex or misbehaved memory allocation.")
      }
    if(C.out$is_error==2)
      stop("in d_dm_dat: misbehaved memory allocation")

    return(C.out$result)
  }




