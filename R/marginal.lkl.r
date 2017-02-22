##' Estimates the marginal likelihood of a model, proceeding by simple Monte-Carlo integration under the prior distribution.
##'
##' The function is a wrapper calling  \code{\link{MCpriorIntFun}} with parameter \code{FUN} set to \code{likelihood}.
##' @title Marginal model likelihood
##' @inheritParams posteriorMCMC
##' @inheritParams MCpriorIntFun
##' @param dat The angular data set relative to which the marginal model likelihood is to be computed 
##' @param likelihood The likelihood function of the model.
##' See \code{\link{posteriorMCMC}} for the required format.
##' @param displ logical. If \code{TRUE}, a plot is produced, showing the temporal evolution of the cumulative mean, with  approximate confidence intervals of \eqn{+/-2}  estimated standard errors.
##' @param precision  the desired relative precision. See
##' \code{\link{MCpriorIntFun}}.
##' @return  The list returned by \code{\link{MCpriorIntFun}}. The estimate is the list's element  named \code{emp.mean}.
##' @note The estimated standard deviations of the estimates produced by this function should be handled with care:For "larger" models than the Pairwise Beta or the NL models,
##' the  likelihood may  have
##' infinite second moment under the prior distribution.  In such a case, 
##' it is recommended to resort to more sophisticated integration methods,
##' \emph{e.g.} by sampling from a mixture of the prior and the
##' posterior distributions. See the reference below for more details.
##' @export
##' @references KASS, R. and  RAFTERY, A. (1995). Bayes factors. \emph{Journal of the american statistical association , 773-795}.
##' @examples
##' \dontrun{
##'   lklNL=  marginal.lkl(dat=Leeds,
##'                  likelihood=dnestlog,
##'                  prior=prior.nl,
##'                  Nsim=20e+3,
##'                  displ=TRUE,
##'                  Hpar=nl.Hpar,
##'                 )
##'}
##'
##' @seealso \code{\link{marginal.lkl.pb}}, \code{\link{marginal.lkl.nl}} for direct  use with  the implemented models.
marginal.lkl <-
function(dat,
         likelihood,
         prior,
         Nsim=300 ,
         displ=TRUE,
         Hpar,
         Nsim.min=Nsim,
         precision=0,
         show.progress = floor(seq(1, Nsim, length.out = 20 ) )
         )
  {
    intern.fun=function(param)
      {
        likelihood(x=dat, par=param,
                   log = FALSE, vectorial = FALSE)
      }
    ## MC integration 
    mc.res=MCpriorIntFun(Nsim=Nsim,
      prior=prior,
      Hpar=Hpar,
      dimData= ncol(dat),
      FUN=intern.fun,
      store=TRUE,
      Nsim.min=Nsim.min, precision=precision,
      show.progress=show.progress
      )
  
    ## end.
    ## check convergence of marginal likelihood estimators
    if(displ )
      {
        nsim = mc.res$nsim
        est.mean = cumsum(mc.res$stored.vals) / (1:nsim)
        est.err = sqrt( cumsum( (mc.res$stored.vals- est.mean)^2) )/
          (1:nsim)

        dev.new()
        par(mfrow=c(1,1))
   
        plot(1:nsim,est.mean,type="l",
             xlab="mc iterations",
             ylab="likelihood",
             main="", lwd=1.5,
             ylim=c(min( est.mean,
               (est.mean-2*est.err)[floor(nsim/2):nsim]),
               max(est.mean, (est.mean+2*est.err)[floor(nsim/2):nsim])
               )
             )
        lines(1:nsim, est.mean + 2* est.err,col=gray(0.5))
        lines(1:nsim , est.mean -  2* est.err,col=gray(0.5))
        ## legend("topright",
        ##         legend=c("mean likelihood","2*stand. deviation of  mean estimator"),
        ##         #"2*stand.dev. of successive points" 
        ##         lwd=c(1,1),
        ##         col=c("black","gold") ,
        ##         cex=0.7)
          
      }
  
    return(mc.res)

  }

