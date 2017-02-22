##' Simple Monte-Carlo sampler approximating the integral of \code{FUN} with respect to the prior distribution. 
##'
##' The algorithm exits after \eqn{n} iterations,
##' based on the following stopping rule :
##'  \eqn{n} is the minimum number of iteration, greater than
##' \code{Nsim.min},  such that the relative
##' error is less than the specified \code{precision}. 
##' \deqn{ max (est.esterr(n)/ |est.mean(n)| ) \le \epsilon ,} where
##' \eqn{est.mean(n)} is the estimated mean of \code{FUN} at time 
##' \eqn{n}, \eqn{est.err(n)} is the estimated  standard
##' deviation of the estimate:
##' \eqn{est.err(n) = \sqrt{est.var(n)/(nsim-1)} }.
##' The empirical variance is computed component-wise and the maximum
##' over the parameters' components is considered.
##'
##' The algorithm exits in any case after \code{Nsim} iterations, if the above condition is not fulfilled before this time. 
##' @title Generic Monte-Carlo integration of a function under the prior distribution
##' @param Nsim Maximum number of iterations
##' @inheritParams posteriorMCMC 
##' @param dimData The dimension of the model's \emph{sample} space,
##' on which the parameter's dimension may depend.
##' Passed to \code{prior} inside \code{MCintegrateFun}
##' @param FUN A function to be integrated. It may return a vector or an array.
##' @param store Should the successive evaluations of \code{FUN} be stored ?
##' @param show.progress same as in \code{\link{posteriorMCMC}}
##' @param Nsim.min The minimum number of iterations to be performed.
##' @param precision The desired relative precision \eqn{\epsilon}.
##' See \bold{Details} below.  
##' @param ... Additional arguments to be passed to \code{FUN}. 
##' @return  A list made of
##' \itemize{
##' \item \code{stored.vals} : A matrix with \code{nsim} rows and
##' \code{length(FUN(par))} columns.
##' \item \code{elapsed} : The time elapsed during the computation.
##' \item \code{nsim} : The number of iterations performed
##' \item \code{emp.mean} : The desired integral estimate: the empirical mean.
##' \item \code{emp.stdev} : The empirical standard deviation of the sample.
##' \item \code{est.error} : The estimated standard deviation of the estimate (\emph{i.e.} \eqn{emp.stdev/\sqrt(nsim)}).
##' \item \code{not.finite} : The number of non-finite values obtained (and discarded) when evaluating \code{FUN(par,...)}
##' }
##' @author Anne Sabourin
##' @export
MCpriorIntFun <-
function(Nsim=200,
         prior,
         Hpar,
         dimData,
         FUN=function(par,...){as.vector(par)},
         store=TRUE,
         show.progress = floor(seq(1, Nsim, length.out = 20 ) ),
         Nsim.min=Nsim,
         precision = 0,
         ...)
  {
 
############ intialize  ############

    start.time=proc.time()
    not.finite=0
    param = prior(type = "r", n=1, Hpar=Hpar, dimData=dimData)
    temp.res=FUN(param,...)
    dim.res=dim(temp.res)
    
    if(is.null(dim.res) || (sum(dim.res!=1) ==1)  )
      {
        emp.mean=rep(0,length(temp.res))
      }
    else
      {
        store=FALSE
        emp.mean=array(0,dim=dim.res)
      }
    emp.variance= emp.mean
    emp.variance.unNorm=emp.variance

    if(store)
      {
        stored.vals=matrix(0,nrow=Nsim,ncol=length(emp.mean))
      }
      
################ start MC ########
    nsim=1
    while((nsim<=Nsim) &&
          ( (nsim<=Nsim.min) ||
           (max( sqrt(emp.variance/(nsim-1)) /
                abs(emp.mean) ) > precision) )
          )
      {
        ## show progression      
        if(any(nsim==show.progress))
          {
            cat(paste((nsim-1), "iterations done", "\n", sep = " " ))
          }


        flag=TRUE
        count = 0
        while(flag & (count<=50))
          {
            param = prior(type = "r", n=1, Hpar=Hpar, dimData=dimData)
            temp.res=FUN(param,...)
            flag = (any(sapply(as.vector(temp.res),
                               function(x){ ! is.finite(x) } ) ) )
            if(flag)
              {
                not.finite = not.finite+1
              }
            count = count+1
          }
        if(flag)
          stop("more than 50 non finite values produced in a row")
        
        cur.res=temp.res 
        
        new.emp.mean=emp.mean+1/nsim*(cur.res-emp.mean)
        
        emp.variance.unNorm=emp.variance.unNorm +
          (cur.res-new.emp.mean)* (cur.res- emp.mean)
        emp.variance = emp.variance.unNorm/(nsim-1)
        emp.mean = new.emp.mean
        
        if(store)
          {
            stored.vals[nsim,]= as.vector(cur.res)
          }
        nsim=nsim+1
      }
######### end MC #########
    end.time = proc.time()
    elapsed=end.time-start.time
    print(elapsed)
    if(store)
      {
        returned.vals=stored.vals[1:(nsim-1),]
      }
    else
      {
        returned.vals=0
      }
    
    return(list( stored.vals= returned.vals,
                    elapsed=elapsed,
                    nsim = nsim-1,
                    emp.mean=emp.mean,
                    emp.stdev=sqrt(emp.variance),
                    est.error=sqrt(emp.variance/(nsim-1)),
                    not.finite = not.finite))
    
  }
