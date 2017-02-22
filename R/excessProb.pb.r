

##' Simple MC integration on the simplex for joint excess probability,
##' in the PB model.
##'
##' @title Estimates the probability of joint excess, given a PB parameter.
##' @param par the DM parameter, as a list
##' @param thres the multivariate threshold
##' @param precision The desired relative precision of the estimate.
##' @param Nmin The number of MC iterations to be performed
##' @param displ logical: should convergence diagnostic plots be issued ?
##' @param add logical: should the plot be added to a current one ? 
##' @return a list made of \describe{
##' \item{mean}{The mean estimate from the MC sample}
##' \item{esterr}{The estimated standard deviation of the estimator}
##' \item{estsd}{The estimated standard deviation of the MC sample}
##' }
##' @export
##' @keywords internal
excessProb.condit.pb <- function(par=c(0.8,1,2,3),
                                  thres=rep(500,5),
                                 precision=0.1, Nmin=200,
                                 displ=FALSE, add=FALSE)
  {
    ## res <- rep(0,N)
    ## ws <- rpairbeta(N,par=par,dimData=3)
    ## res <- apply(ws, 1, function(w){3*min(w/thres)})
    i <- 0
    cond <- FALSE
    mean <- 0
    Mi <- 0
    si <- 0
    res <- double(2*Nmin)
    while((i<Nmin) || !cond)
      {
        i <- i+1
        if(i>length(res))
          res <- c(res,double(i))

        repeat{    
          w <- as.vector(rpairbeta(n=1,par=par, dimData=3))
          if(all(w>0)) break
        }
        xi <- 3*min(w/thres)
        res[i] <- xi
        delta <- xi-mean
        mean <- mean + delta/i
        Mi <- Mi + delta*(xi-mean)
        si <- Mi/i
        cond <-  (si/(i*mean^2) < precision^2 || i> 1e+5)
      }
    N <- i
    res <- res[1:N]
    if(displ)
      {
        cummean <- cumsum(res)/(1:N)
        estsd <- sqrt(cumsum((res-cummean)^2)/(1:N) )
        esterr <- estsd/sqrt(1:N)

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
    ##    return(list(mean=cummean[N], esterr=esterr[N], estsd=estsd[N]))
    return(c(mean, sqrt(si/N) , sqrt(si)))
  }



## tt <- excessProb.condit.pb(par=c(0.8,1,2,3),
##                                   thres=rep(500,3),
##                                  precision=0.01, Nmin=200,
##                                  displ=T, add=FALSE)

#NULL 

##' Double Monte-Carlo integration. 
##'
##' @title Estimates the probability of joint excess (Frechet margins)
##' @inheritParams posteriorMean
##' @param Nmin.intern The minimum number of MC iteration in the internal loop (excess probability, conditional to a parameter).
##' @param precision The desired precision for the internal MC estimate
##' @param post.sample The  posterior sample.
##' @param thres A multivariate threshold
##' @param known.par Logical
##' @param true.par The true parameter from which the data are issued.
##' @return A list made of \describe{
##' \item{whole}{ A vector of estimated excess probabilities, one for each element of the thinned posterior sample.}
##' \item{mean}{the estimated threshold excess probability: mean estimate.}
##' \item{esterr}{The estimated standard deviation of the mean estimate
##' (where the Monte-Carlo error is neglected)}
##' \item{estsd}{The estimated standard deviation of the posterior sample (where the Monte-Carlo error is neglected)}
##' \item{lowquants}{The three lower \eqn{0.1} quantiles of, respectively, the conditional mean estimates and of the upper and lower bounds of the Gaussian (centered) \eqn{80} \% confidence intervals around the conditional estimates. }
##' \item{upquants}{The three upper \eqn{0.9} quantiles}
##' \item{true.est}{the mean estimate conditional to the true parameter:
##' a vector of size three: the mean estimate , and the latter +/- the standard deviation of the estimate}
##' }
##' @export
excessProb.pb <- function( post.sample,
                          Nmin.intern=100, precision=0.05,
                          from=NULL,to=NULL, thin=100,
                          displ=FALSE,
                          thres=rep(500,5), known.par= FALSE,
                          true.par)
  {
    reslist <- posteriorMean(post.sample=post.sample,
                             from=from,to=to, thin=thin,
                             FUN=excessProb.condit.pb,
                             Nmin=Nmin.intern, precision=precision,
                             thres=thres,
                             displ=FALSE)

    res <- reslist$values
    N <- ncol(res)
    cummean <- cumsum(res[1,])/(1:N) 
    estsd <- sqrt(cumsum((res[1,]-cummean)^2)/(1:N)) 
    esterr <- estsd/sqrt(1:N)
    ymax= max(res[1,]) ## max(cummean+1.1*estsd)
    ymin=min(res[1,])  ##min(cummean-1.1*estsd)
    if(displ)
      {
        plot(1:N, cummean, ylim=range(res[1,]), type="l", lwd=2) ##mean estimate

        polygon(c(1:N, N:1),
                c(cummean+qnorm(0.9)*estsd,
                  rev(cummean-qnorm(0.9)*estsd)), col=gray(0.8))
        lines(1:N, cummean,   lwd=2) ##mean estimate
      }
    if(known.par)
      {
        true.est <- excessProb.condit.pb ( par=true.par,Nmin=Nmin.intern,
                                          precision=precision,
                                          thres=thres,
                                          displ=FALSE)
        if(displ){    
          polygon(c(1,N,N,1),
                  c(true.est[1]+ qnorm(0.9)*true.est[2],
                    true.est[1]+ qnorm(0.9)*true.est[2],
                    true.est[1] - qnorm(0.9)*true.est[2],
                    true.est[1] - qnorm(0.9)*true.est[2]),
                  density=10, col="red")
        
          abline(h=true.est[1], col="red", lwd=2 )
        }
      }
       
    else{ true.est <- NULL }
    sorted <- ##sort.int(res[,1], decreasing=TRUE, index.return=TRUE)
      sort(res[1,], decreasing=TRUE)
    sortedUp <- sort(res[1,] + qnorm(0.9)* res[2,], decreasing=TRUE)
    sortedLow <- sort(res[1,] - qnorm(0.9)* res[2,], decreasing=TRUE)
    
    upquant <- sorted[ceiling(10/100*N)] ##sorted$x[ceiling(10/100*N)]
    upquantUp <- sortedUp[ceiling(10/100*N)]
    upquantLow <- sortedLow[ceiling(10/100*N)]
    
    ## polygon(c(1,N,N,1), c(upquantUp,upquantUp,upquantLow,upquantLow),
    ##         col="blue", density=10)
    if(displ){
      abline(h=upquant, col="blue",  lwd=2)
    }
    lowquant <- sorted[floor(90/100*N)]
    lowquantUp <- sortedUp[floor(90/100*N)]
    lowquantLow <- sortedLow[floor(90/100*N)]

    ## polygon(c(1,N,N,1), c(lowquantUp,lowquantUp,lowquantLow,lowquantLow),
    ##         col="blue", density=10)
    if(displ){
      abline(h=lowquant, col="blue", lwd=2)
      lines(1:N, cummean,   lwd=2) ##mean estimate 
         
      if(known.par)
        {
          legend("topright",
                 legend=c("true", "posterior mean",
                   "posterior 0.1/0.9 quantiles",
                   "posterior 0.1/0.9 Gaussian quantiles" ),
                 lwd=c(2,2,2,3),
                 col=c("red", "black", "blue", gray(0.5))
                 )
        }
      else
        {
          legend("topright", legend=c( "posterior mean",
                               "posterior 0.1/0.9 quantiles",
                               "posterior 0.1/0.9 Gaussian quantiles" ),
                 lwd=c(2,2,4),
                 col=c( "black", "blue", gray(0.5))
                 )
      
        }
    }
    return(list(whole=res, mean=cummean[N], esterr=esterr[N],
                estsd=estsd[N],
                lowquants=c( lowquant,lowquantLow,lowquantUp),
                upquants=c(upquant,upquantLow,upquantUp),
                true.est=true.est))
  }

     
