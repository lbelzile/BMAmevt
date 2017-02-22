##' @title Maximum likelihood optimization
##' @param data The angular data to be used for inference
##' @param model A list made of \describe{
##'\item{likelihood}{The likelihood function, see \code{\link{dpairbeta}}
##' for a template}.
##' \item{npar}{The length of the parameter vector}
##' }
##' @param init NULL or a real vector of size \code{model$npar} giving the initial values for \code{link{par}}.
##' @param maxit maximum number of iterations to be performed by
##' function \code{optim}
##' @param method The method to be used by \code{optim}
##' @param hess logical: should an approximation of the hessian be performed ?
##' @param link the link function from the natural marginal parameter spaces to the real line.
##' @param unlink the inverse link function. If \code{x} is any real number, then \code{unlink(x)} should be in the admissible range for the likelihood function and the prior function.
##' @return The list returned by \code{optim} and the AIC and BIC criteria
##' @export
maxLikelihood <-
  function(data, model,## = list(likelihood, npar), 
           init = NULL, maxit = 500, method="L-BFGS-B", hess = T, link, unlink)
  ## @param prior The prior function (see \code{\link{prior.pb}} for a template) for generating initial parameters in case the initial value results in non finite log-likelihood.
## @param Hpar the prior hyper parameters

{
  p = dim(data)[2]
  ndat=dim(data)[1]
  npar <- model$npar
  
  lhood <- function(vec)
    {
      -model$likelihood(x=data, par=unlink(vec),
                        log=TRUE, vectorial=FALSE)
    }


  if(is.null(init))
    init <- rep(0,model$npar)
  count <- 0
  converged <- FALSE
  while(!converged & count<20 )
    {
      count <- count+1
      opt.test <- tryCatch(optim(init, lhood, method = method,
                                 control =list(maxit = maxit,trace=0 ),
                                 hessian = hess),
                           error=function(e){
                           #  print(e)
                             return(list(convergence=100))
                           })
      if(as.integer(opt.test$convergence)>0){
        init <- rnorm(model$npar,mean=0,sd=1)
      }
       
      else{
        converged <- TRUE
        opt <- opt.test
      }

    }
  rtn <- list()
  if(as.integer(opt.test$convergence)>0){
    rtn$message <- "optimisation failed"
    rtn$counts <- count
    rtn$convergence <- 100
    rtn$linkedpar <- rep(0,npar)
    rtn$par <- unlink(rep(0,npar))
    rtn$value <- -Inf
    rtn$aic <- rtn$aicc <- rtn$bic <- Inf
    rtn$linkedHessian <- 1
    return(rtn)

  }
  
  rtn$message <-  opt$message
  rtn$counts <-  opt$counts 
  rtn$convergence <- opt$convergence
  rtn$linkedpar <- opt$par
  rtn$par <- unlink(opt$par)
  rtn$value <- opt$value
  rtn$aic <- 2* (opt$value + npar)
  rtn$aicc <- 2* (opt$value + npar) +2*npar*(npar+1)/(ndat-npar-1 )
  rtn$bic <-  2*opt$value + npar*log(ndat)

  if(!hess)
    {
      rtn$linkedHessian <- 1
    }
  else
    {
      rtn$linkedHessian <- opt$hessian
      tryCatch(expr={
      asympt.variance <- chol2inv(chol(opt$hessian))
      rtn$asympt.variance <- asympt.variance
      rtn$linked.esterr <- sqrt(diag(asympt.variance))},
      error=function(e){## rtn$asympt.variance <- NULL
                        ## rtn$linked.esterr <- NULL
                      return(rtn)} )
    }
        
  return(rtn)
}
##?optim
