##' The exponent function \eqn{V} for a max-stable variable \eqn{M} is such that \eqn{P(M<x) = exp(-V(x))}
##'
##' @title Exponent function in the NL  model.
##' @param par The parameter for the NL  distribution,
##' respectively of length two or four.
##' @param x A vector of three extended positive real numbers
##' @return the value of \eqn{V(x)} for \eqn{x=thres}.  
##' @export
expfunction.nl <- function(par=c(0.3,0.4,0.5,0.6),
                        x=10*rep(1,3))
  {
    alpha=par[1]
    beta=par[-1]

    U1 = ( x[1]^(-1./(alpha*beta[1]) ) + 
      x[2]^(-1./(alpha*beta[1]) ) )^(beta[1])
    
    U2 = ( x[1]^(-1./(alpha*beta[2]) ) + 
      x[3]^(-1./(alpha*beta[2]) ) )^(beta[2])

    
    U3 = ( x[2]^(-1./(alpha*beta[3]) ) + 
      x[3]^(-1./(alpha*beta[3]) ) )^(beta[3])

    return( 2^(-alpha) * ( (U1 + U2 + U3)^(alpha)  ) )
    
  }


##' @title Probability of joint threshold excess in the NL model
##' @param par The Nested logistic parameter: of length four.
##' @param thres a positive vector of size three.
##' @return The approximate probability of joint excess, valid when at least one coordinate of \code{thres} is large
##' @export
excessProb.condit.nl <- function(par=c(0.3,0.4,0.5,0.6),
                                  thres=rep(100,3))
  {
    expfun <- expfunction.nl

    zeros <- which(thres==0)
    if(length(zeros)==3)
          stop(" x must have at least one positive element")
    if(length(zeros)==2)
      {
        x <- thres
        x[zeros] <- Inf
        return(expfun(par=par, x=x))
      }
    
    if(length(zeros)==1)
      {
        nonzeros <- c(1,2,3)[-zeros]
        x <- thres
        x[zeros] <- Inf
        T0 <- expfun(par=par,x=x)
          
        x1 <- x
        x1[nonzeros[1]] <- Inf
        T1 <- expfun(par=par, x=x1)
        
        x2 <- x
        x2[nonzeros[2]] <- Inf
        T2 <- expfun(par=par, x=x2)

        return( T1 + T2 - T0 )
      }
    
    T1 <- expfun(par=par,x=thres)
    
    T2 <- expfun(par=par,x=c(thres[1], Inf,Inf)) +
      expfun(par=par,x=c(Inf, thres[2],Inf)) +
        expfun(par=par,x=c(Inf, Inf,  thres[3]))

    T3 <- expfun(par=par,x=c(thres[1], thres[2],Inf)) +
      expfun(par=par,x=c(thres[1], Inf, thres[3])) +
        expfun(par=par,x=c(Inf, thres[2],  thres[3]))
    
    return( T1 + T2 - T3 )
    
  }
                                
##' @title Posterior distribution the probability of joint threshold excess, in the NL model.
##' @param post.sample The posterior sample, as returned by \code{posteriorMCMC}
##' @inheritParams posteriorMean
##' @inheritParams excessProb.condit.nl
##' @param known.par logical. Is the true parameter known ?
##' @param true.par The true parameter, only used  if \code{known.par=TRUE} 
##' @return A list made of \describe{
##' \item{whole}{The output of \code{posteriorMean} called with \code{FUN=excessProb.condit.nl}.}
##' \item{mean}{The posterior mean of the excess probability}
##' \item{esterr}{The standard deviation of the mean estimator}
##' \item{estsd}{The standard deviation of the excess probability,
##' in the posterior sample. }
##' \item{lowquant}{The lower 0.1 quantile of the empirical posterior distribution of the excess probability }
##' \item{upquant}{The upper 0.1 quantile of the empirical posterior distribution of the excess probability }
##' \item{true}{\code{NULL} if \code{known.par=FALSE}, otherwise the excess probability in the true model.}
##' }
##' @export
excessProb.nl <- function(post.sample,
##                          model="trinestlog",
                          from=NULL,to=NULL, thin=100,
                          thres=rep(100,3),
                          known.par= FALSE,
                            true.par,
                          displ=FALSE)
  {
    reslist <- posteriorMean(post.sample=post.sample,
                             from=from,to=to, thin=thin,
                             FUN=excessProb.condit.nl,
                             displ=FALSE,
                             thres=thres
                             )
                                        #browser()    
    res <- as.vector(reslist$values)
    N <- length(res)
    cummean <- cumsum(res)/(1:N) 
    estsd <- sqrt(cumsum((res-cummean)^2)/(1:N)) 
    esterr <- estsd/sqrt(1:N)
    ymax= max(res) ## max(cummean+1.1*estsd)
    ymin=min(res)  ##min(cummean-1.1*estsd)

    if(displ){
    
      plot(1:N, cummean, ylim=range(res), type="l", lwd=2) ##mean estimate

      polygon(c(1:N, N:1),
              c(cummean+qnorm(0.9)*estsd,
                rev(cummean-qnorm(0.9)*estsd)), col=gray(0.8))
      lines(1:N, cummean,   lwd=2) ##mean estimate
    }
    if(known.par)
      {

        true <- excessProb.condit.nl ( par=true.par,
                                      thres=thres
                                     )
        if(displ)
          abline(h=true, col="red", lwd=2 )
      }
    else
      true <- NULL
    
    sorted <- ##sort.int(res[,1], decreasing=TRUE, index.return=TRUE)
      sort(res, decreasing=TRUE)
    upquant <- sorted[ceiling(10/100*N)]
    lowquant <- sorted[floor(90/100*N)]
    if(displ){
      abline(h=lowquant, col="blue", lwd=2)
      abline(h=upquant, col="blue", lwd=2)
      if(known.par)
        {
          legend("topright", legend=c("true", "posterior mean",
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
                lowquant= lowquant,
                upquant=upquant,
                true=true))


  }
