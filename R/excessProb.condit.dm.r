
##' simple MC integration on the simplex.
##'
##' @title Probability of joint threshold exceedance, in the Dirichlet Mixture model,  given a DM parameter.
##' @param N The number of MC iterations to be performed
##' @param par the DM parameter, as a list
##' @param thres the multivariate threshold
##' @param plot logical: should convergence diagnostic plots be issued ?
##' @param add logical: should the plot be added to a current one ? 
##' @return a list made of \describe{
##' \item{mean}{the mean estimate from the MC sample}
##' \item{esterr}{the estimated standard deviation of the estimator}
##' \item{estsd}{The estimated standard deviation of the MC sample}
##' }
##' @export
excessProb.condit.dm <- function(N=100, par=get("dm.expar.D3k3"),
                                 thres=rep(100,3),
                                 plot=FALSE, add=FALSE)
  {
    res <- rep(0,N)
    for(i in 1:N)
      {
        w <- as.vector(rdirimix(n=1,par=par))
        res[i] <-nrow(par$Mu )*min(w/thres)
      }
    cummean <- cumsum(res)/(1:N)
    estsd <- sqrt(cumsum((res-cummean)^2)/(1:N) )
    esterr <- estsd/sqrt(1:N)
    if(plot)
      {
        ymax=max(cummean+1.1*estsd)
        ymin=min(cummean-1.1*estsd)
        if(!add)
          plot(1:N, cummean, ylim=c(ymin,ymax), type="l",col="blue")
        else
          lines(1:N, cummean,col="blue")
        
        lines(1:N, cummean+esterr, col="blue", lty=2)
        lines(1:N, cummean-esterr,col="blue", lty=2)
        ## lines(1:N, cummean+estsd, col="red")
        ## lines(1:N, cummean-estsd,col="red")
      }
    return(list(mean=cummean[N], esterr=esterr[N], estsd=estsd[N]))
  }

