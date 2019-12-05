##' Approximation of a model marginal likelihood by Laplace method.
##'
##' The posterior mode is either supplied, or approximated by numerical optimization. For an introduction about Laplace's method, see \emph{e.g.}
##' Kass and Raftery, 1995 and the references therein.
##' @title Laplace approximation of  a model marginal likelihood by Laplace approximation. 
##' @param mode The parameter vector (on the \dQuote{unlinked} scale, \emph{i.e.} before transformation to the real line)
##' which maximizes the posterior density, or \code{NULL}. 
##' @param npar The size of the parameter vector. Default to four. 
##' @param likelihood The likelihood function, \emph{e.g.} \code{\link{dpairbeta}} or \code{\link{dnestlog}}
##' @param prior The prior density (takes an \dQuote{unlinked} parameter as argument and returns the density of the \code{linked} parameter)
##' @param Hpar The prior hyper parameter list.
##' @param data The angular dataset
##' @param link The link function, from the \dQuote{classical} or \dQuote{unlinked} parametrization onto the real line. (\emph{e.g.} \code{log} for the PB model, an \code{logit} for the NL model) 
##' @param unlink The inverse link function (\emph{e.g.} \code{exp} for the PB model and \code{invlogit} for the NL model)
##' @param method The optimization method to be used. Default to \code{"L-BFGS-B"}.
##' @return A list made of \describe{
##'\item{mode}{the parameter (on the unlinked scale) deemed to maximize the posterior density. This is equal to the argument if the latter is not null.}
##' \item{value}{The value of the posterior, evaluated at \code{mode}.}
##' \item{laplace.llh}{The logarithm of the estimated marginal likelihood}
##' \item{invHess}{The inverse of the estimated hessian matrix at \code{mode}}
##' }
##' @references KASS, R.E. and RAFTERY, A.E. (1995). Bayes Factors.
##' \emph{Journal of the American Statistical Association,
##' Vol. 90, No.430}
##' @export
laplace.evt <- function(mode=NULL, npar=4, likelihood, prior, Hpar,
                        data, link, unlink, method="L-BFGS-B")

  {
    if(!is.null(mode))
      npar <- length(mode)

    post.fun <- function(par)
      {
        likelihood(x=data, par=unlink(par),  log=TRUE, vectorial=FALSE)+
          prior(type="d",par=unlink(par),log=TRUE, Hpar=Hpar)
      }

    if(is.null(mode)){
        count <- 0
        converged <- FALSE
        init <- rep(0,npar)
        while(!converged & count<20 ){
            count <- count+1
            opt.test <-tryCatch( optim(par=init, fn=post.fun, gr = NULL,
                                       method=method,
                                       control = list(fnscale=-1),
                                       hessian =  FALSE),
                                error=function(e) return(
                                  list(convergence=100))
                                )
                                
            if(as.integer(opt.test$convergence)>0)
              init <- rnorm(npar,mean=0,sd=1)
            
            else{
              converged <- TRUE
              opt <- opt.test
            }

          }
      if(as.integer(opt.test$convergence)){
        res <- list(mode=unlink(rep(0,npar)), value=-Inf,
                    laplace.llh = -Inf,
                    invHess = 1)
        return(res)
      }     
        linkedmode <- opt$par
        mode <- unlink(opt$par)
      }
    
    else
      linkedmode <- link(mode)
    
    fmode <- post.fun(linkedmode)
    
    hessian.mode <- tryCatch(optimHess(par=linkedmode,
                                       fn=post.fun, gr = NULL),
                             error=function(e) return( 1))
    hessian.mode <- 0.5*(hessian.mode+t(hessian.mode))
    
    ## Inverse.test <- try(VarCov.t <-  -solve(hessian.mode), ##-chol2inv(chol(hessian.mode)),#, 
    ##                     silent = TRUE)
    Inverse.test <- try(invHess.t <-  chol2inv(chol(-hessian.mode)),#, 
                        silent = TRUE)

    if (!inherits(Inverse.test, "try-error")) {
      invHess <- Inverse.test
      diag(invHess) <- ifelse(diag(invHess) <= 0, .Machine$double.eps, 
                             diag(invHess))
        }
    else {
      cat("\nWARNING: singular Hessian matrix in Laplace.evt\n ")
      cat("\nNB: Identity matrix is used instead of the inverse ,
 consider the result with care.\n")
      invHess <- diag(npar)
    }
    laplace.llh <- NA
#    options(warn = -1)
    laplace.test <- try(LML0 <- npar/2 * log(2 * pi) + 0.5 * 
                    log(det(invHess)) +
                    as.vector(fmode), 
                    silent = TRUE)
    if (is.finite(laplace.test[1])) 
      laplace.llh <- laplace.test[1]
 #   options(warn = 0)
    res <- list(mode=mode, value=fmode,
                    laplace.llh = laplace.llh,
                    invHess = invHess,
                hessian=hessian.mode)
  return(res)
  }

